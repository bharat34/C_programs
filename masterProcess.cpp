#include "masterProcess.hpp"
# This program inititlaize all the objects for the initial pdb input file
# such as charnge, van der waal radii, connectivity information
# residues information, secondary structure information, persisistance length
# and gyration, dihedral information 
masterProcess::masterProcess(){
	time_t t1,t2;
	time(&t1);

	//I want to have that number of TOP_STRUCTS_MAX which be >=100 and can be equally divided in 
	//the processors available, so in Boltzmann MC, each processor has equal no. of structs.
	TOP_STRUCTS_MAX=100;

	//Read the input protein struct.
	char file_name[20]="job_input.pdb";
	inputStruct=pdb(file_name);
	ressz=inputStruct.residueSize();	//no of residues.
	res=new residue[ressz]; assert(res!=NULL);

	//Constucting the residue information.
	for(int i=0; i<inputStruct.residueSize(); i++){
		res[i].create(inputStruct, i+1);
		res[i].setPolar(inputStruct);
	}

	//Setting the charge, Rst & epsilon for all atoms from the AminoRad.
	AminoRad arad;
	for(int i=1; i<=inputStruct.size(); i++){
		arad.setParameter(inputStruct[i].getResidueSequence(), ressz, inputStruct[i].getAtomName().c_str(),\
			inputStruct[i].getResidueName().c_str());
		inputStruct[i].Ep=arad.epsilon;
		inputStruct[i].R=arad.r;
		inputStruct[i].charge=arad.charge;
		inputStruct[i].atomType = arad.atomType;
	}

	//Joining adjacent residues at a time to get the bonded information.
	unsigned pdbSize= inputStruct.size();
        unsigned startAtom=0;
	for(int i=0; i<inputStruct.residueSize()-1; i++){
		res[i].join(res[i+1]);
		startAtom= (i<ressz-2) ? res[i+2].getStartAtom() : 0;
		res[i].getBondedAtomList(bondList,angleList,dihedralList,impDihedralList,nonBondList,pdbSize,startAtom);
		res[i].setRamachandranAngles(res[i+1]);
		res[i].remove(res[i+1]);
	}
		res[inputStruct.residueSize()-1].getBondedAtomList(bondList,angleList,dihedralList,impDihedralList,nonBondList,pdbSize,startAtom);// this statement was found missing, inserted on 13-03-06

	char file_name2[20]="dat/parm94.dat";
	Parameter parm(file_name2);
	setBondParameter(inputStruct, parm, bondList);
	setAngleParameter(inputStruct, parm, angleList);
	setDihedralParameter(inputStruct, parm, dihedralList);
	setIDihedralParameter(inputStruct, parm, impDihedralList);
	setNonBondParameter(inputStruct, nonBondList); //dielectric provided during calculation of electrostatic E.

	//Reading the secondary structure information.
	FILE *ssfp = fopen("job_input.ss","r"); assert(ssfp!=NULL);
	totalSecStrs=0; //total number of secondary structs.
        char line[80];
        while(fgets(line,80,ssfp)!=NULL){
	        sscanf(line,"%c%d%d",&secStrInfo[totalSecStrs].type, &secStrInfo[totalSecStrs].start, \
				&secStrInfo[totalSecStrs].end);
                totalSecStrs++;
        }
        fclose(ssfp);

	//Changing the secondary struct limits if a loop is of less than 4 residues.
	for(int i=0; i<(totalSecStrs-1); ++i){
		//if(secStrInfo[i+1].start<=(secStrInfo[i].end+1)){
		if(secStrInfo[i+1].start<=(secStrInfo[i].end)){   ///modified for taking care of continuous helix - sheet
			cout<<"No AA in loop region: "<<secStrInfo[i].end+1<<" "<<secStrInfo[i+1].start-1<<endl;
			delete this;
			exit(-1);
		}

		int length=secStrInfo[i+1].start-secStrInfo[i].end-1;
		bool alternate=true;
		while(length<4){
			if(alternate){
				--secStrInfo[i].end;
				++length;	
				alternate=false;
			}
			else{
				++secStrInfo[i+1].start;
				++length;	
				alternate=true;
			}
		}
	}

	//finding the max length of secondary structure.
	maxSecLength=0;
	for(int i=0; i<totalSecStrs; ++i){
		int length=secStrInfo[i].end-secStrInfo[i].start+1;
		cout<<i<<"th Secondary Struct=> start:"<<secStrInfo[i].start<<" end:"<<secStrInfo[i].end<<endl;
		if(maxSecLength<length) maxSecLength=length;
	}
	cout<<"maxSecLength: "<<maxSecLength<<endl;
	
	//Assigning the sec str info to residues.
        for(int i=0; i<totalSecStrs; ++i){
	        for(unsigned j=secStrInfo[i].start; j<=secStrInfo[i].end; j++){
	                res[j-1].flagSS=secStrInfo[i].type;
			res[j-1].nSS=i;
                }
        }

	//Generating the rotatable dihedrals indices(numbering is 1,2,3.....)
	totalRotateDD=0;
	int diheds[7];
	for(int i=0 ; i<(totalSecStrs-1); ++i){  	//(totalSecStrs-1) loops.
		for(int j=0; j<7; ++j) diheds[j]=0;

		//find the rotatable diheds in for the ith loop, given [i+1]secStr start & [i]secStr end.
		rotatable_diheds(secStrInfo[i+1].start, secStrInfo[i].end, diheds);

		//storing & printing the rotatable dihedrals.
		int t=0;
		cout<<"Rotatable diheds for "<<(i+1)<<"th loop are: "<<endl;
		while(diheds[t] && t<7){
			if(diheds[t]%2) cout<<"\t\t"<<(diheds[t]+1)/2<<"-->phi"<<endl;	//2n-1
			else cout<<"\t\t"<<diheds[t]/2<<"-->psi"<<endl;			//2n

			rotateDD[totalRotateDD]=diheds[t];
			++t;
			++totalRotateDD;
		}	
	}

	//Setting AA_index hash_map.
	AA_index["ALA"]=ALA; AA_index["ARG"]=ARG; AA_index["ASN"]=ASN;
	AA_index["ASP"]=ASP; AA_index["CYS"]=CYS; AA_index["GLN"]=GLN;
	AA_index["GLU"]=GLU; AA_index["GLY"]=GLY; AA_index["HIS"]=HIS;
	AA_index["ILE"]=ILE; AA_index["LEU"]=LEU; AA_index["LYS"]=LYS;
	AA_index["MET"]=MET; AA_index["PHE"]=PHE; AA_index["PRO"]=PRO;
	AA_index["SER"]=SER; AA_index["THR"]=THR; AA_index["TRP"]=TRP;
	AA_index["TYR"]=TYR; AA_index["VAL"]=VAL;

	nextStructSeqNo=0; postLoopCount=0; postPLCount=0; postRGCount=0; 
	dihedralIndices=new int[totalRotateDD];
	dihedralValues=new int[totalRotateDD];
	assert(dihedralIndices!=NULL); assert(dihedralValues!=NULL);

	maxTS=(long)pow(2.0, totalRotateDD);//==> 2^#of rotatable dihedrals..

	time(&t2);
	cout<<"Master process boot time: "<<difftime(t2,t1)<<endl;
}

//Find the seven or less dihedrals for rotation.
void masterProcess::rotatable_diheds(const unsigned start2, const unsigned end1, int *diheds){

	int length=start2-end1-1; 		///using the residue amino acid
	int resIndex[MAX_LOOP_LENGTH];		///store the index of the dihedral; index=0(nonproline) and 1(proline)
	int countDihed=0;			///total no of dihed available for rotation

	for(int i=0;i<length;i++){
		if(res[end1+i].getResidueName()=="PRO"){
			countDihed++;
			resIndex[i]=1;
		}
		else{
			resIndex[i]=0;
			countDihed+=2;
		}
	}
	int count,endRotate1,startRotate2; ///endRotate1 contains 3 
	if(countDihed<=7){
		count=0;
		for(int i=0;i<length;i++){
			if(resIndex[i]==0){	///not a proline
				diheds[count]=(2*(end1+1+i))-1;////2n - 1
				count++;
				diheds[count]=(2*(end1+1+i));///2n
				count++;
			}
			else{ 			///for a proline
				diheds[count]=(2*(end1+1+i)); /////2n for proline
				count++;
			}
		}
	}
	else if(countDihed>7){
		count=0;
		for(int i=0; i<length; i++){  	///select first 3 dihed
			endRotate1=end1+1+i;
			if(resIndex[i]==0){			///not a proline
				if(count==3) break;
				diheds[count]=(2*(end1+1+i))-1;	///2n - 1, n is amino acid
				count++;

				if(count==3) break;
				diheds[count]=(2*(end1+1+i));	///2n
				count++;
			}
			else{ 					///for a proline
				if(count==3) break;
				diheds[count]=(2*(end1+1+i)); 	///2n for proline
				count++;
			}
		}
		count=6;
		for(int i=length-1;i>=0;i--){		///select last 3 dihed
			startRotate2=end1+1+i;
			if(resIndex[i]==0){			///not a proline
				if(count==3) break;
				diheds[count]=(2*(end1+1+i));	///2n  psi first(reverse)
				count--;

				if(count==3) break;
				diheds[count]=(2*(end1+1+i)-1);	///2n -1
				count--;
			}
			else{ 					///for a proline
				if(count==3) break;
				diheds[count]=(2*(end1+1+i)); 	///2n for proline ie consider psi
				count--;
			}
		}
		/// Changed as phi of proline was selected for the case of 1ybz(18804848)
		if(  ((diheds[2]+diheds[4]) % 2)==0  ) 
			diheds[3]=(diheds[2]+diheds[4])/2;
		else
			diheds[3]=(diheds[2]+diheds[4])/2 + 1;

		if((diheds[3]==diheds[2])||(diheds[3]==diheds[4])){
			delete this;
			exit(-1);
			cout<<"Problem in rotatable AAs"<<endl;
			delete this;
			exit(-1);
		}	
			
		//if((startRotate2 - endRotate1) > 1) diheds[3]=2*(endRotate1 + (startRotate2-endRotate1)/2);
		//else diheds[3]=2*(startRotate2);	
	}
}

double masterProcess::calPerLen(const pdb &currentTS){
	double maxLen=0, len=0, lenTotMax=0;
	for(int i=0;i<ressz-1;i++){
		unsigned flag=0;
		for(int j=i+1;j<ressz;j++){
			len=currentTS[res[i].getN()].getCoordinates().distance(currentTS[res[j].getN()].getCoordinates());
			if(flag!=0){
				if(len<maxLen) break;
				else{
					maxLen=len;
					if(maxLen>lenTotMax) lenTotMax=maxLen;
				}
			}
			if(flag==0){
				maxLen=len;
				if(maxLen>lenTotMax) lenTotMax=maxLen;
				flag=1;
			}
		}
	}

/*
	for(unsigned i=0; i<totalSecStrs; i++){
		for(unsigned j=i; j<totalSecStrs; j++){
			len=(currentTS[res[secStrInfo[i].start-1].getCA()].pt - currentTS[res[secStrInfo[j].end-1].\
						getCA()].pt).mod_square();
			if(len>maxLen) maxLen=len;
			else{
				if(maxLen > lenTotMax) lenTotMax=maxLen;
				maxLen=0;
				break;
			}
		}
	}
	return(sqrt(lenTotMax));
*/

	return(lenTotMax);
}

//Radius of Gyration Filter.
double masterProcess::gyrationRadius(const pdb &currentTS){

	double totAtWt=0, atWt, totAtWtX=0, totAtWtY=0, totAtWtZ=0;
	double Cx, Cy, Cz; ////center of Gravity 
	double atomRG, totRG=0; 	
	unsigned i;
	for(i=1; i<=currentTS.size(); i++){

		atWt = currentTS[i].getAtomWt();
		totAtWt += atWt;

		totAtWtX = totAtWtX + ( atWt * currentTS[i].getCoordinates().x); ///Wt * x coordinates
		totAtWtY = totAtWtY + ( atWt * currentTS[i].getCoordinates().y); ///Wt * y coordinates
		totAtWtZ = totAtWtZ + ( atWt * currentTS[i].getCoordinates().z); ///Wt * z coordinates
	}
	Cx = totAtWtX / totAtWt;		
	Cy = totAtWtY / totAtWt;		
	Cz = totAtWtZ / totAtWt;		
		
	for(i=1; i<=currentTS.size(); i++){
		atomRG =( ((currentTS[i].getCoordinates().x - Cx) * (currentTS[i].getCoordinates().x - Cx)) + \
			((currentTS[i].getCoordinates().y - Cy) * (currentTS[i].getCoordinates().y - Cy)) + \
			((currentTS[i].getCoordinates().z - Cz) * (currentTS[i].getCoordinates().z - Cz)) );
		totRG += atomRG;
	}

	return(sqrt(totRG/(currentTS.size()+1)));
}


//If possiblity of the generation of a new trial struct then the generated trial struct is 
//returned(referenced) in newTS else -1 is returned.
int masterProcess::genStruct(pdb &newTS){
	bool newTSFound=false;
	bool postLoopFilter=false;

	//nextStructSeqNo is the seq no of the current structure to be generated.
	while(!newTSFound && (nextStructSeqNo<maxTS)){
		unsigned temp=nextStructSeqNo;
		//Getting the indices of dihedrals required.
		for(int i=0; i<totalRotateDD; ++i){
			dihedralIndices[i]=temp%2;
			temp/=2;
		}

		postLoopFilter=true;
		++nextStructSeqNo;

		for(int i=0; postLoopFilter&&(i<totalRotateDD); ++i){
			if(rotateDD[i]%2){ //PHI.
				int aa=(rotateDD[i]+1)/2;
				int allowedDihedI=AA_index[res[aa-1].getResidueName().c_str()]; 
								//for the residue resName, allowedDihedI 
								//gives the allowed phi & psi
								//values from maps AMINO_PHI/PSI_DIHEDS.

				int phi=AMINO_PHI_DIHEDS[allowedDihedI][dihedralIndices[i]];

				if(phi == -1) postLoopFilter=false;
				else{
					//Have a valid phi value.
					dihedralValues[i]=phi;
				}	
			}
			else{ //PSI.
				int aa=rotateDD[i]/2;
				int allowedDihedI=AA_index[res[aa-1].getResidueName().c_str()]; 
								//for the residue resName, allowedDihedI 
								//gives the allowed phi & psi
								//values from maps AMINO_PHI/PSI_DIHEDS.

				int psi=AMINO_PSI_DIHEDS[allowedDihedI][dihedralIndices[i]];

				if(psi == -1) postLoopFilter=false;
				else{
					//Have a valid psi value.
					dihedralValues[i]=psi;
				}	
			}
		}	

		if(postLoopFilter){		//We have a valid set of dihedral values for the loops.
			++postLoopCount;

			//set the dihedral values for the trial struct.
			for(int i=0; i<totalRotateDD; ++i){
				if(rotateDD[i]%2){ //PHI.
					int resNo=(rotateDD[i]+1)/2 -1;
					res[resNo].setPhiAngle(newTS, dihedralValues[i]);
				}	
				else{	//PSI
					int resNo=rotateDD[i]/2 -1 ;
					res[resNo].setPsiAngle(newTS, dihedralValues[i]);
				}	
			}

/*
			//Write the trial structs.
			char name[20];
			sprintf(name, "%dTS.pdb", nextStructSeqNo-1);
			newTS.write(name);
*/
			if(getTotalSecStrs()==2)   ////
				return 0;

			double perLength=calPerLen(newTS);
			//perlength upper limit is #of residues in longest sec struct * 1.5 *4.
			if((perLength>=20.2)&&(perLength<=(maxSecLength*6.0))){
				++postPLCount;

				//if longest sec struct <=20 only then use RG filter.
				if(maxSecLength<=20){
					double RG=gyrationRadius(newTS);
					double Y1=0.3952*pow(ressz,0.6)+4.2572;
					double Y2=0.3952*pow(ressz,0.6)+15.2572;
					if((Y1<=RG)&&(RG<=Y2)){
						++postRGCount;
						newTSFound=true;

/*
						//Write the candidate trial structs.
						char name[20];
						sprintf(name, "%drg.pdb", nextStructSeqNo-1);
						newTS.write(name);
*/
					}
				}
				else{//RG filter is not to be applied.
					--postRGCount;
					newTSFound=true;
				}
			}
		}
	}

	if(newTSFound) return 0; 
	else return -1;				//no more trial structs can be generated.
}

/* A slaveProcess reports a trial structure(post clash removal & energy minimization) */
void masterProcess::postClashMini(const double *coordsRecv){
	const double energy=coordsRecv[0];
	const double *coords=&(coordsRecv[1]);
	const int nAtoms=inputStruct.size();

	list<pdbDescription>::iterator it=topStructsList.begin();
	int precedingCount=0;
	//list is in the ascending order of energy.
	while((it!=topStructsList.end()) && (it->energy<energy)){
		++it;
		++precedingCount;
	}
	//cout<<"precedingCount: "<<precedingCount<<endl;

	if(precedingCount<TOP_STRUCTS_MAX){	//the new structure should be in top structs.
		pdbDescription pd;
		pd.energy=energy;
		pd.coordinates=new double[3*nAtoms];
		assert(pd.coordinates != NULL);

		for(int i=0; i<3*nAtoms; ++i) pd.coordinates[i]=coords[i];
		topStructsList.insert(it, pd);
	}
	
	if(topStructsList.size()>TOP_STRUCTS_MAX){	//the list has more than TOP_STRUCTS_MAX, remove last.
		pdbDescription &pd=topStructsList.back();
		delete[] pd.coordinates;
		topStructsList.pop_back();
	}
}

void masterProcess::writeTopStructs(){
	char name[10];
	list<pdbDescription>::iterator it=topStructsList.begin();
	const int nAtoms=inputStruct.size();

	for(int i=0; (i<TOP_STRUCTS_MAX)&&(it!=topStructsList.end()); ++i,++it){
		for(int j=0; j<nAtoms; ++j){
			inputStruct[j+1].pt.x=it->coordinates[3*j];
			inputStruct[j+1].pt.y=it->coordinates[3*j+1];
			inputStruct[j+1].pt.z=it->coordinates[3*j+2];
		}
		sprintf(name, "%dP.pdb", i);
		inputStruct.write(name);
	}
}

masterProcess::~masterProcess(){
	cout<<"Total trial structs generated: "<<nextStructSeqNo<<endl;
	cout<<"Number of structures post Loop filter: "<<postLoopCount<<endl;
	cout<<"Number of structures post persistence length filter: "<<postPLCount<<endl;
	cout<<"Number of structures post radius of gyration filter: "<<postRGCount<<endl;

	delete[] res;
	delete[] dihedralIndices;
	delete[] dihedralValues;

	list<pdbDescription>::iterator it=topStructsList.begin();
	while(it!=topStructsList.end()){
		delete[] it->coordinates;
		++it;
	}

	topStructsList.clear();
	bondList.clear();
	angleList.clear();
	dihedralList.clear();
	impDihedralList.clear();
	nonBondList.clear();
}
