#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <set>
#include<sys/stat.h>
#include "SecStr.hpp"
#include "makefragment.cpp"
#include "stch.cpp"
# Last Modified by Bharat Lakhani on (01-23-13)
# program reads all the small patches of the pdb file to 
# make the full fragment. This a parrelized MPI code and not
# limited to specific number of processors.
using namespace std;

#define WORKTAG 1
#define RESULTTAG 2
#define DIETAG 3
#define MAX_LEN 100
#define MIN_LEN 20


void master();
void slave(int rank);

typedef struct{
	char pdb[5];
	char Id[2];
        int start;
        int end;
	string dirname;
	vector < vector <string> > fragmentsFile1;
	vector < vector <string> > fragmentsFile2;
}AlignmentInfo;


int file_exist (char *filename)
{
  struct stat   buffer;
  return (stat (filename, &buffer) == 0);
}



char filename[MIN_LEN];
char jobID[MIN_LEN];


int main(int argc,char **argv){
        MPI::Init(argc, argv);
        int myRank=MPI::COMM_WORLD.Get_rank();

        if(!myRank) cout<<"Running bhageerath-H Secondary Protocol: "<<endl;

        time_t t1, t2;

        time(&t1);

	strcpy(filename,argv[1]);
        strcpy(jobID,argv[2]);

        if(!myRank) master();
        else slave(myRank);

        time(&t2);

        if(!myRank) cout<<"Time taken in bhageerath-H Secondary Protocol: "<<difftime(t2,t1)<<endl;

        MPI::Finalize();
        return 0;
                              }

void master() {
 
	int nTasks=MPI::COMM_WORLD.Get_size();
        MPI::Status status;
        MPI::Request req;

	FILE *f1=fopen(filename,"r");
	list<AlignmentInfo> Alignment;  
	list<AlignmentInfo>::iterator it1,it2,itnew;
	list<AlignmentInfo> AlignmentDummy;  
	it1=Alignment.begin();

	char store[1000]; int start,end; char pdb[5]; char Id[2], dirname[15];
	
	FILE *writefiles;
	int writefiles_c = 0;
	
	int mynode=1, length;

	while(fgets(store,1000,f1)!=NULL) 
	{
						
			sscanf(store," %s  %s  %d  %d  %*d  %s",pdb,Id,&start,&end,dirname);
			
			sprintf(store,"%s/%s/job_input.aa.B99990001.pdb",jobID,dirname);
			if (file_exist(store) != 1) continue;

			sprintf(store,"%s/%s/job_input.aa.B99990002.pdb",jobID,dirname);
			if (file_exist(store) != 1) continue;

		        it1=Alignment.begin();

	while(it1!=Alignment.end() && ((it1->start < start && it1->end < end ) || (it1->start > start && it1->end < end) || (it1->start == start && it1->end < end)|| (it1->start < start && it1->end == end )|| (it1->start == start && it1->end == end )))
			 ++it1;

			 sprintf(store,"%s/%s/job_input.aa.B99990001.pdb",jobID,dirname);
			 char deletelines[strlen(store)+126];
			 sprintf(deletelines,"sed -i '/^END/d;/^TER/d' %s/%s/job_input.aa.B99990001.pdb",jobID,dirname);
			 system(deletelines);
			 vector <string> tmpFragments;
			 fstream fp_pdb;
			 AlignmentInfo align;
                         int FindTheLimit = 1;
                         int swap =0;
			 if (file_exist(store) == 1)
			 {
				fp_pdb.open(store);
				string str;
				getline (fp_pdb,str);
				getline (fp_pdb,str);
				getline (fp_pdb,str);
				while (getline (fp_pdb,str))
				{
					istringstream buffer(str.substr(22,10));
					 buffer >> swap;
					if (FindTheLimit!=swap )
					{
						align.fragmentsFile1.push_back(tmpFragments);
						tmpFragments.clear();
						FindTheLimit=swap;
						str.insert(str.length(),"\n");
						tmpFragments.push_back(str);
					
					}else
					 {
						str.insert(str.length(),"\n");
						tmpFragments.push_back(str);
					 }

				}
						align.fragmentsFile1.push_back(tmpFragments);
				
			 }
			 tmpFragments.clear();
			 fp_pdb.close();
			 sprintf(store,"%s/%s/job_input.aa.B99990002.pdb",jobID,dirname);
			 sprintf(deletelines,"sed -i '/^END/d;/^TER/d' %s/%s/job_input.aa.B99990002.pdb",jobID,dirname);
			 system(deletelines);
			 if (file_exist(store) == 1)
			 {
				fstream fp_pdb;
				fp_pdb.open(store);
				string str;
				getline (fp_pdb,str);
				getline (fp_pdb,str);
				getline (fp_pdb,str);
				int FindTheLimit = 1;
				int swap =0;
				while (getline (fp_pdb,str))
				{
					istringstream buffer(str.substr(22,10));
					 buffer >> swap;
					if (FindTheLimit!=swap )
					{
						align.fragmentsFile2.push_back(tmpFragments);		
						tmpFragments.clear();
						FindTheLimit=swap;
						str.insert(str.length(),"\n");
						tmpFragments.push_back(str);
					
					}else
					 {
						str.insert(str.length(),"\n");
						tmpFragments.push_back(str);
					 }
				}
						align.fragmentsFile2.push_back(tmpFragments);		
				
			 }
			 tmpFragments.clear();
			 align.start=start;
			 align.end=end;
			 strcpy(align.pdb,pdb);
			 strcpy(align.Id,Id);
			 align.dirname=dirname;
			 Alignment.insert(it1, align);	
	 }
/*
		
		it1=Alignment.begin();
		while(it1 != Alignment.end())
         {
		cout<<it1->start<<"	"<<it1->end<<endl;
		it1++;
	}
*/
		 it1=Alignment.begin();
			int count=0;

	 while(it1 != Alignment.end())
         {

                count=0;
                itnew=Alignment.begin();
                while (itnew !=Alignment.end())
                {
                        if(it1->start==itnew->start && it1->end == itnew->end && count <=5)
                        {

                                count++;
                                cout<<it1->start<<"     "<<it1->end<<"   "<<"  "<<itnew->start<<"  "<<itnew->end<<endl;
                        }
                      else if(it1->start==itnew->start && it1->end == itnew->end && count >5) {itnew=Alignment.erase(itnew);   }
                        itnew++;

                }


                it1++;
        }
	itnew=Alignment.begin();
	while(itnew!=Alignment.end()){ cout<<"*****"<<itnew->start<<"  "<<itnew->end<<endl; itnew++;}

		
	//Make the Combination store the position in the list
			int combination = 0;
			it1=Alignment.begin();
			signed temp1, temp2;

	 while(it1 != Alignment.end())     
	 {
			continue;
			 cout<<"Find the patch for : "<<it1->start<<"   "<<it1->end<<endl;
	 		 AlignmentDummy.assign(Alignment.begin(),Alignment.end());
			 it2 = AlignmentDummy.begin();
		         temp1 = it1->start; temp2 = it1->end;
			 int position=0;
			 list<unsigned int> positionlist;
			 list<unsigned int>::iterator pos1;
	 while(it2 != AlignmentDummy.end())
	 {

	 if(temp1 == it2->start && temp2 == it2->end) 
	 { 
		it2++; 
		position++;
		continue; 
	 }
	 if(((temp1 > it2->start && temp2 == it2->end)||(temp1 > it2->start && temp2 < it2->end) || (temp1 == it2->start && temp2 < it2->end) || (temp1 < it2->start && temp2 < it2->end && temp2 > it2->start)) &&((it2->end-temp2 >=10) || (temp1-it2->start >=10))) 
	{
	    
		if(temp1 > it2->start) temp1 = it2->start; 
		if(temp2 < it2->end) temp2 = it2->end;
		positionlist.push_back(position);
	//	cout<<position<<"	"<<it2->start<<"	"<<it2->end<<endl;
	}
		it2++;
		position++;
	}
	
	pos1=positionlist.begin();

	vector< vector<string> >::iterator it_fragmentsFiles;
	vector<string>::iterator it_profile;
	vector <string> tmpFragments1;
	vector <string> tmpFragments2;
	vector <string> dummyFragments1;
	vector <string> dummyFragments2;

	for (it_fragmentsFiles = (*it1).fragmentsFile1.begin(); it_fragmentsFiles!=(*it1).fragmentsFile1.end(); it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) tmpFragments1.push_back(*it_profile);

	for (it_fragmentsFiles = (*it1).fragmentsFile2.begin(); it_fragmentsFiles!=(*it1).fragmentsFile2.end(); it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) tmpFragments2.push_back(*it_profile);


	//Make Combination with each position and send it to slave node
	while(positionlist.size()!=0)
	{ 

	combination++;
	if(positionlist.size() > 3 ) {  positionlist.pop_front(); continue; }

	temp1 = it1->start; temp2 = it1->end;
	string patchinfo;
	char yyy[500];
	sprintf(yyy,"Patched For	:%d	%d	%s 	%s 	%s\n",temp1,temp2,it1->pdb,it1->Id,it1->dirname);
	cout<<"Patched For	:"<<temp1<<"	"<<temp2<<" "<<it1->pdb<<" "<<it1->Id<<"	"<<it1->dirname<<endl;
	patchinfo = yyy;

	dummyFragments1.assign(tmpFragments1.begin(),tmpFragments1.end());
	dummyFragments2.assign(tmpFragments2.begin(),tmpFragments2.end());
	int howmanypatches = 0;
	
	for(pos1=positionlist.begin(); pos1!=positionlist.end(); pos1++) 
	{

	it2=Alignment.begin();
	advance(it2,*pos1);

	if(((temp1 > it2->start && temp2 == it2->end)||(temp1 > it2->start && temp2 < it2->end) || (temp1 == it2->start && temp2 < it2->end) || (temp1 < it2->start && temp2 < it2->end && temp2 > it2->start)) &&((it2->end-temp2 >=10) || (temp1-it2->start >=10)))
	{
	//cout<<it2->start<<"		"<<it2->end<<"	"<<it2->dirname<<endl;
	howmanypatches++;
	

	if(temp1 > it2->start && (temp1-it2->start >=10)) 
	{
			 

	vector <string> templist;
		
	//for (it_fragmentsFiles = (*it2).fragmentsFile1.begin() ; it_fragmentsFiles!= (*it2).fragmentsFile1.begin() + ( temp1 -  it2->start +1); it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) templist.push_back(*it_profile);
	 for (it_fragmentsFiles = (*it2).fragmentsFile1.begin() ; it_fragmentsFiles != ( (*it2).fragmentsFile1.begin() +  (temp1 - it2->start )); it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) templist.push_back(*it_profile);
	templist.push_back("TER\n");
	it_profile = dummyFragments1.begin();
	dummyFragments1.insert(it_profile,templist.begin(),templist.end());
	templist.clear();

	//for (it_fragmentsFiles = (*it2).fragmentsFile2.begin() ; it_fragmentsFiles!= (*it2).fragmentsFile2.begin() + ( temp1 - it2->start + 1); it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) templist.push_back(*it_profile);
	for (it_fragmentsFiles = (*it2).fragmentsFile2.begin() ; it_fragmentsFiles !=(*it2).fragmentsFile2.begin() + ( temp1 - it2->start ); it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) templist.push_back(*it_profile);
	templist.push_back("TER\n");
	it_profile = dummyFragments2.begin();
	dummyFragments2.insert(it_profile,templist.begin(),templist.end());
	templist.clear();
	temp1 = it2->start; 
	}

        if(temp2 < it2->end && (it2->end-temp2 >=10)) 
	{ 

	dummyFragments1.push_back("TER\n");
	//for (it_fragmentsFiles = (*it2).fragmentsFile1.begin() + ( temp2 - it2->start + 1 ); it_fragmentsFiles!= ((*it2).fragmentsFile1.begin() + (it2->end-it2->start + 1)); it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) dummyFragments1.push_back(*it_profile);
	for (it_fragmentsFiles = (*it2).fragmentsFile1.begin() + ( temp2 -it2->start +1   ); it_fragmentsFiles != (*it2).fragmentsFile1.end(); it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) dummyFragments1.push_back(*it_profile);
	dummyFragments2.push_back("TER\n");
	
	//for (it_fragmentsFiles = (*it2).fragmentsFile2.begin() + ( temp2 - it2->start + 1 ); it_fragmentsFiles!= ((*it2).fragmentsFile2.begin() + (it2->end-it2->start + 1)); it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) dummyFragments2.push_back(*it_profile);
	for (it_fragmentsFiles = (*it2).fragmentsFile2.begin() + ( temp2 -it2->start +1  ); it_fragmentsFiles != (*it2).fragmentsFile2.end() ; it_fragmentsFiles++) for (it_profile = (*it_fragmentsFiles).begin(); it_profile != (*it_fragmentsFiles).end(); it_profile++ ) dummyFragments2.push_back(*it_profile);

	temp2 = it2->end; 
	}

        sprintf(yyy,"New Patch	:	%d	%d	%s	%s 	%s \n",temp1,temp2,it2->pdb,it2->Id,it2->dirname);
        cout<<"New Patch	:"<<temp1<<"	"<<temp2<<" "<<it2->pdb<<" "<<it2->Id<<"	"<<it2->dirname<<endl;
	patchinfo.append(yyy);	
	}
	}

				string sendpdb1;
				string sendpdb2;
for(int kk=0; kk < dummyFragments1.size(); kk++)  sendpdb1.append(dummyFragments1[kk],0,dummyFragments1[kk].length());
for(int kk=0; kk < dummyFragments2.size(); kk++)  sendpdb2.append(dummyFragments2[kk],0,dummyFragments2[kk].length());
		
		
	 if(howmanypatches != 0)
	 {
		continue;
	 if(mynode != nTasks) 
	 {
		length = sendpdb1.size() + 1;
		req = MPI::COMM_WORLD.Isend(&length, 1, MPI::INT, mynode, WORKTAG); req.Wait(status);
		req = MPI::COMM_WORLD.Isend(sendpdb1.c_str(),length, MPI::CHAR, mynode, WORKTAG); req.Wait(status);
		length = patchinfo.size() + 1;
		req = MPI::COMM_WORLD.Isend(&length, 1, MPI::INT, mynode, WORKTAG); req.Wait(status);
		req = MPI::COMM_WORLD.Isend(patchinfo.c_str(),length, MPI::CHAR, mynode, WORKTAG); req.Wait(status);
		mynode++;

	 if(mynode != nTasks) 
	 {
		length = sendpdb2.size() + 1;
		req = MPI::COMM_WORLD.Isend(&length, 1, MPI::INT, mynode, WORKTAG); req.Wait(status);
		req = MPI::COMM_WORLD.Isend(sendpdb2.c_str(),length, MPI::CHAR, mynode, WORKTAG); req.Wait(status);
		length = patchinfo.size() + 1;
		req = MPI::COMM_WORLD.Isend(&length, 1, MPI::INT, mynode, WORKTAG); req.Wait(status);
		req = MPI::COMM_WORLD.Isend(patchinfo.c_str(),length, MPI::CHAR, mynode, WORKTAG); req.Wait(status);
		mynode++;	
	 }
	 else 
	 {
	
		char name[40];
		sprintf(name,"%dT.pdb",writefiles_c);
		FILE *writefiles = fopen(name,"w");
	

	 	int recvlentgh1;
		MPI::COMM_WORLD.Recv(&recvlentgh1, 1, MPI::INT, MPI::ANY_SOURCE, RESULTTAG, status);
		char getpdb1[recvlentgh1];
		getpdb1[recvlentgh1] = '\0';
	 	MPI::COMM_WORLD.Recv(&getpdb1, recvlentgh1, MPI::CHAR, status.Get_source(), RESULTTAG, status);
/*	
		fprintf(writefiles,"%s",getpdb1);
		fclose(writefiles);
*/	
		
			
		length = sendpdb2.size() + 1;
		MPI::COMM_WORLD.Send(&length, 1, MPI::INT, status.Get_source(), WORKTAG);
		MPI::COMM_WORLD.Send(sendpdb2.c_str(), length, MPI_CHAR,status.Get_source(), WORKTAG);

		length = patchinfo.size() + 1;
		MPI::COMM_WORLD.Send(&length, 1, MPI::INT, status.Get_source(), WORKTAG);
		MPI::COMM_WORLD.Send(patchinfo.c_str(), length, MPI_CHAR,status.Get_source(), WORKTAG);
		writefiles_c++;
	 }
	 }				
	else 
	{
	
		char name1[40];
                sprintf(name1,"%dT.pdb",writefiles_c);
                FILE *writefiles = fopen(name1,"w");
	

		int recvlentgh1;
		MPI::COMM_WORLD.Recv(&recvlentgh1, 1, MPI::INT, MPI::ANY_SOURCE, RESULTTAG, status);
		char getpdb1[recvlentgh1];
		getpdb1[recvlentgh1] = '\0';
		MPI::COMM_WORLD.Recv(&getpdb1, recvlentgh1, MPI::CHAR, status.Get_source(), RESULTTAG, status);
/*	
		fprintf(writefiles,"%s",getpdb1);
                fclose(writefiles);
	*/	writefiles_c++;
	

		length = sendpdb1.size() + 1;
		MPI::COMM_WORLD.Send(&length, 1, MPI::INT, status.Get_source(), WORKTAG);
		MPI::COMM_WORLD.Send(sendpdb1.c_str(), length, MPI_CHAR,status.Get_source(), WORKTAG);
		length = patchinfo.size() + 1;
		MPI::COMM_WORLD.Send(&length, 1, MPI::INT, status.Get_source(), WORKTAG);
		MPI::COMM_WORLD.Send(patchinfo.c_str(), length, MPI_CHAR,status.Get_source(), WORKTAG);

	
		char name2[40];
                sprintf(name2,"%dT.pdb",writefiles_c);
                writefiles = fopen(name2,"w");
	

		int recvlentgh2;
		MPI::COMM_WORLD.Recv(&recvlentgh2, 1, MPI::INT, MPI::ANY_SOURCE, RESULTTAG, status);
		char getpdb2[recvlentgh2];
		getpdb2[recvlentgh2] = '\0';
		MPI::COMM_WORLD.Recv(&getpdb2, recvlentgh2, MPI::CHAR, status.Get_source(), RESULTTAG, status);
	/*
		fprintf(writefiles,"%s",getpdb2);
                fclose(writefiles);
	*/	writefiles_c++;
	

		length = sendpdb2.size();
		MPI::COMM_WORLD.Send(&length, 1, MPI::INT, status.Get_source(), WORKTAG);
		MPI::COMM_WORLD.Send(sendpdb2.c_str(), length, MPI_CHAR,status.Get_source(), WORKTAG);
		length = patchinfo.size() + 1;
		MPI::COMM_WORLD.Send(&length, 1, MPI::INT, status.Get_source(), WORKTAG);
		MPI::COMM_WORLD.Send(patchinfo.c_str(), length, MPI_CHAR,status.Get_source(), WORKTAG);

	}
	}
					dummyFragments1.clear();
					dummyFragments2.clear();
					sendpdb1.clear();
					sendpdb2.clear();
					positionlist.pop_front();
	}
		AlignmentDummy.clear();
		positionlist.clear();
		it1++;				
	
	}				
	while(mynode != nTasks) {
                                  length = strlen(filename);
                                  MPI::COMM_WORLD.Send(&length, 1, MPI::INT, mynode, DIETAG);
				  mynode++;
				}
	 mynode=1;
	while(mynode != nTasks) {	
	
		continue;
		char name2[40];
                sprintf(name2,"%dT.pdb",writefiles_c);
                writefiles = fopen(name2,"w");
	
	
		int recvlentgh=0;
		MPI::COMM_WORLD.Recv(&recvlentgh, 1, MPI::INT, MPI::ANY_SOURCE, RESULTTAG, status);
		char getpdb[recvlentgh];
		getpdb[recvlentgh] = '\0';
		MPI::COMM_WORLD.Recv(getpdb, recvlentgh, MPI::CHAR, MPI::ANY_SOURCE, RESULTTAG, status);
	/*
		fprintf(writefiles,"%s",getpdb);
                fclose(writefiles);
	*/	writefiles_c++;
	

		MPI::COMM_WORLD.Send(&length, 1, MPI::INT, status.Get_source(), DIETAG);
		mynode++;
				}
	cout<<"Total Combination made :"<<combination<<endl;
}

void slave(int myRank) 
{	

		MPI::Request reqRecv;
                MPI::Request reqSend;
                MPI::Status status;
                int length=0,len=0;
                int processed_str=0;

		SecStr secStrInfo[MAX_LEN];   

		FILE *ssfp = fopen("job_input.ss","r"); assert(ssfp!=NULL);
 		int totalSecStrs=0; //total number of secondary structs.
        	char line[80];
        	while(fgets(line,80,ssfp)!=NULL){
               	sscanf(line,"%c%d%d",&secStrInfo[totalSecStrs].type, &secStrInfo[totalSecStrs].start, \
                                &secStrInfo[totalSecStrs].end);
                totalSecStrs++;
        					 }
        	fclose(ssfp);

		
		FILE *aafp = fopen("job_input.aa","r"); assert(aafp!=NULL);
        	char sequence[5000];
        	//(fgets(sequence,5000,ssfp));
        	(fgets(sequence,5000,aafp));
        	fclose(aafp);

		 char dirName[80];
                strcpy(dirName,"/home/bhageerathH/PdbInputFiles/"); //pdbinputfiles
		
		string stpdb=stch(dirName,sequence);
		//char stchain[stpdb.size()];
		//strcpy(stchain,stpdb.c_str());
		//int stsize=stpdb.size();


	
		while(true){
			

		reqRecv = MPI::COMM_WORLD.Irecv(&length, 1, MPI::INT, 0, MPI::ANY_TAG);
                reqRecv.Wait(status);


		if(status.Get_tag()==DIETAG)  break;

		++processed_str;

		char getpdb[length];
		getpdb[length]='\0';


                reqRecv = MPI::COMM_WORLD.Irecv(&getpdb, length, MPI::CHAR, 0, MPI::ANY_TAG);
                reqRecv.Wait(status);

		reqRecv = MPI::COMM_WORLD.Irecv(&length, 1, MPI::INT, 0, MPI::ANY_TAG);
                reqRecv.Wait(status);

		char patchinfo[length];
               	patchinfo[length]='\0';

                reqRecv = MPI::COMM_WORLD.Irecv(&patchinfo, length, MPI::CHAR, 0, MPI::ANY_TAG);
                reqRecv.Wait(status);
	
		length = strlen(getpdb);
		char stchain[stpdb.size()];
		strcpy(stchain,stpdb.c_str());
		int stsize=stpdb.size();
		//cout<<stchain<<endl;
		string newpdb = makefragment(totalSecStrs,secStrInfo,getpdb,length,sequence,patchinfo,stchain,stsize);

		string appendpdb;
		
		appendpdb = patchinfo;
		appendpdb.append(newpdb);
		newpdb.clear();

		length = appendpdb.size() + 1;
	

		reqSend=MPI::COMM_WORLD.Isend(&length, 1, MPI::INT, 0, RESULTTAG);
                reqSend.Wait(status);

		
	 	reqSend=MPI::COMM_WORLD.Isend(appendpdb.c_str(), length, MPI::CHAR, 0, RESULTTAG);
                reqSend.Wait(status);

			   }

		cout<<"Structures processed by "<<myRank<<"th slave are: "<<processed_str<<endl;
}
