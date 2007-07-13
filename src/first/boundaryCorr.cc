/*
 Brian Patenaude, parts taken from bet2
 FMRIB Image Analysis Group
 */


//#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "math.h"

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "meshclass/meshclass.h"
#include "shapeModel/shapeModel.h"
using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace mesh;
using namespace shapemodel;


string title="overlapMeasures (Version 1.0) University of Oxford (Brian Patenaude)";
string examples="overlapMeasure [options] -i segImage -l goldStandard -k output.txt ";


Option<bool> verbose(string("-v,--verbose"), false, 
					 string("switch on diagnostic messages"), 
					 false, no_argument);
Option<bool> help(string("-h,--help"), false,
				  string("display this message"),
				  false, no_argument);
Option<bool> inputprob(string("--inputprob"), false,
					   string("Input image is probability - do not do any thresholding or smoothing"),
					   false, no_argument);
Option<bool> outputzmap(string("--zmap"), false,
					   string("outputs zmap"),
					   false, no_argument);

Option<string> inname(string("-i,--in"), string(""),
					  string("filename of input image/mesh/bvars"),
					  true, requires_argument);
Option<string> cmaname(string("-g,--in"), string(""),
					  string("filename of gold standard image"),
					  false, requires_argument);

Option<string> refname(string("-r,--in"), string(""),
					  string("filename of reference image "),
					  false, requires_argument);

Option<int> meshLabel(string("-l,--meshlabel"), 0,
					   string("If loading a vtk mesh, specify label used."),
					   false, requires_argument);
Option<float> thresh(string("-p,--thrsh"), 1.64,
					   string("threshhold for clean up."),
					   false, requires_argument);

Option<string> outname(string("-k,--out"), string(""),
					   string("filename of output mesh"),
					   true, requires_argument);


int nonoptarg;

////////////////////////////////////////////////////////////////////////////
//global variables
float mode(vector<float> vdists,int *maxcount){
	int N=static_cast<int>(vdists.size());
	*maxcount=0;
	float maxlowint=0;
	float bins=256;
	bins=128;
	float binwidth=(vdists.at(N-1)-0)/bins;
//	cout<<"binwidth"<< binwidth<<" "<<vdists.at(N-1)<<" "<<vdists.at(0)<<endl;
		int count=0;	
		int bincount=1;
			float lowint=0;
	for (int i=0;i<N;i++){
	
	
		if (vdists.at(i)<lowint+binwidth){
			count++;
			if (count>(*maxcount)){
				*maxcount=count;
				maxlowint=lowint;
			}
		}else{
			//cout<<count<<" "<<lowint+bincount*binwidth<<" "<<lowint<<endl;
			count=0;
			i--;
			bincount++;
			lowint=bincount*binwidth;
		}
	
	}
	
	//cout<<"bincount "<<bincount<<endl;
	return (maxlowint+binwidth/2.0);
}	
float fullwidthhalfmax(vector<float> vdists,float halfmaxval, float *halfmin,float *halfmax){
	int N=static_cast<int>(vdists.size());
	int maxcount=0;
	float maxlowint=0;
	float bins=256;
	bins=128;
	float binwidth=(vdists.at(N-1)-0)/bins;
	int countprev=0;
	
	bool foundmin=false;
	bool foundmax=false;
	
//	cout<<"binwidth"<< binwidth<<" "<<vdists.at(N-1)<<" "<<vdists.at(0)<<endl;
		int count=0;	
		int bincount=1;
			float lowint=0;
	for (int i=0;i<N;i++){
	
	
		if (vdists.at(i)<lowint+binwidth){
			count++;
			if (count>(maxcount)){
				maxcount=count;
				maxlowint=lowint;
			}
		}else{
			//cout<<count<<" "<<lowint+bincount*binwidth<<" "<<lowint<<endl;
			countprev=count;
			if ((count>halfmaxval)&&(!foundmin)){
				*halfmin=lowint+binwidth/2.0;
				foundmin=true;
			}
			if ((count>halfmaxval)){
				*halfmax=lowint+binwidth/2.0;
			}
			count=0;
			i--;
			bincount++;
			lowint=bincount*binwidth;
		}
	
	}
	
	//cout<<"bincount "<<bincount<<endl;
	return (maxlowint+binwidth/2.0);
}	


void do_work(){
//	if (validate.value()){
		//load gold standard image
		volume<short> mask;
		volume<float> ref;
		read_volume(ref,refname.value());
		read_volume(mask,inname.value());
		int label=meshLabel.value();
	label=10000;
		
		volume<short> gold;
		if (thresh.value()<0){
			read_volume(gold,cmaname.value());
		}
		
		//build intensity vector 
		vector<float> vgraylevels;
		vector <float>::iterator Iter;
		float dist=10000;
		//cout<<
		for (int i=0;i<ref.xsize();i++){
			for (int j=0;j<ref.ysize();j++){
				for (int k=0;k<ref.zsize();k++){
					if ((mask.value(i,j,k)>0)){
						if (mask.value(i,j,k)<label){
							label=mask.value(i,j,k);
						}
						//	cout<<"label is "<<label<<endl;
					}
		}}}
		
		
		for (int i=0;i<ref.xsize();i++){
			for (int j=0;j<ref.ysize();j++){
				for (int k=0;k<ref.zsize();k++){
				
					if (mask.value(i,j,k)==label){
						dist=ref.value(i,j,k);
					
						if (vgraylevels.empty()){
							//	cout<<"empty"<<endl;
							vgraylevels.push_back(dist);
						}else if (dist>=vgraylevels.back()){
							//	cout<<"back"<<endl;
							vgraylevels.push_back(dist);
						}else {
							//	cout<<"insert"<<endl;
							for ( Iter = vgraylevels.begin( ) ; Iter !=vgraylevels.end( ) ; Iter++ ){
								if (dist<*Iter){
									
									vgraylevels.insert(Iter,dist);
									break;
								}
								
							}
							
							
						}
						
					}
				}
			}
		}
		int maxcount;
		float modeval=mode(vgraylevels,&maxcount);
		//cout<<"modeval "<<modeval<<" "<<maxcount<<endl;
		
		//now find full width half maximum
		float halfmaxval=maxcount/2.0;
	
		
		float halfmin,halfmax;
		
		
		float modeval2=fullwidthhalfmax(vgraylevels,halfmaxval,&halfmin, &halfmax);
	
		//cout<<"modeval "<<modeval2<<" "<<halfmin<<" "<<halfmax<<" "<<maxcount<<endl;
		float mean=(halfmin+halfmax)/2;
		float sdev=abs(halfmax-halfmin)/2.35;
		//cout<<"sdev "<<sdev<<endl;
		volume<short> maskcorr;
		volume<float> maskz;
		maskz=ref*0;
		maskcorr=mask;
		
		
		
		
		float zthresh=thresh.value();
		
		for (int i=0;i<ref.xsize();i++){
			for (int j=0;j<ref.ysize();j++){
				for (int k=0;k<ref.zsize();k++){
					float z=0.0;
					if (mask.value(i,j,k)==(label+100  )){
						z=(ref.value(i,j,k)-mean)/sdev;
						maskz.value(i,j,k)=z;
						if (zthresh>=0){
							if (abs(z)>zthresh){
								maskcorr.value(i,j,k)=0;
							}else{
								maskcorr.value(i,j,k)=label;
							}
						}else{
							//this is for boundary corrected dice
							maskcorr.value(i,j,k);
						}
					}
				}
			}
		}
		
		
		
			//this is used for moorph
		cout<<label<<endl;
		
	
		save_volume(maskcorr,outname.value()+"_lcorr");
		if (outputzmap.value()){
			save_volume(maskz,outname.value()+"_zcorr");
		}
}
		
int main(int argc,char *argv[])
{
	
	Tracer tr("main");
	OptionParser options(title, examples);
	
	try {
		// must include all wanted options here (the order determines how
		//  the help message is printed)
		options.add(inname);
		options.add(cmaname);
		options.add(refname);


		options.add(outname);
		options.add(outputzmap);
			options.add(thresh);
			options.add(meshLabel);
	
		
		options.add(inputprob);
		options.add(verbose);
		
		options.add(help);
		
		nonoptarg = options.parse_command_line(argc, argv);
		
		// line below stops the program if the help was requested or 
		//  a compulsory option was not set
		if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
		{
			options.usage();
			exit(EXIT_FAILURE);
		}
		
		// Call the local functions
		do_work();
	}  catch(X_OptionError& e) {
		options.usage();
		cerr << endl << e.what() << endl;
		exit(EXIT_FAILURE);
	} catch(std::exception &e) {
		cerr << e.what() << endl;
	} 
	
	return 0;// do_work(argc,argv);
}

