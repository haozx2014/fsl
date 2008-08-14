/*
 *  shapeModel.cpp
 *  
 *
 *  Created by Brian Patenaude on 23/06/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "shapeModel.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
//#include <vector>

#include <cmath>
#include <algorithm>
#include "math.h"
using namespace std;
namespace SHAPE_MODEL_NAME{
	
	shapeModel::shapeModel()
	{
	
	}
	shapeModel::shapeModel( const vector<float> & mshape, const vector< vector<float> > & modesshape, const vector<float> & se, \
							const vector<float> & ishape, const vector< vector<float> > & modesint, const vector<float> & ie, \
							const int & M, const vector<float> & errs)
{
		smean=mshape;
		smodes=modesshape;
		seigs=se;
		sqrtseigs=se;
		for (vector<float>::iterator i=sqrtseigs.begin();i!=sqrtseigs.end();i++)
			*i = sqrt(*i);
		imean=ishape;
		imodes=modesint;
		ieigs=ie;
		NumberOfSubjects=M;
		Errs=errs;
		USE_COND=false;
		MODE_FOUND=false;
		mode=0;
		
		
}

shapeModel::shapeModel( const vector<float> & mshape, const vector< vector<float> > & modesshape, const vector<float> & se, \
						const vector<float> & ishape, const vector< vector<float> > & modesint, const vector< vector<float> > & Iprec, const vector<float> & ie,\
						const int & M, const vector<float> & errs, const vector< vector<unsigned int> > & cellsin, const vector<int> & vlabels )
{
	labels=vlabels;
	smean=mshape;
	smodes=modesshape;
	seigs=se;
	sqrtseigs=se;
	for (vector<float>::iterator i=sqrtseigs.begin();i!=sqrtseigs.end();i++)
		*i = sqrt(*i);
	imean=ishape;
	imodes=modesint;
	i_precision=Iprec;
	ieigs=ie;
	NumberOfSubjects=M;
	Errs=errs;
	cells=cellsin;
	USE_COND=false;
	MODE_FOUND=false;
	mode=0;
	
	//keeps tracks of neighbourin triangles
	for (unsigned int i=0; i<static_cast<unsigned int>(smean.size()/3); i++)
	{
		vector<unsigned int> inds;
		int count=0;
		for ( vector< vector<unsigned int> >::iterator j=cells.begin(); j!=cells.end();j++,count++)
		{
			for ( vector<unsigned int>::iterator k= j->begin(); k!=j->end();k++)
			{
				if ((*k)==i)
				{
					inds.push_back(count); 
					break;
				}
			}
		}
		localTri.push_back(inds);
	}
	
	
}

shapeModel::shapeModel( const vector<float> & mshape, const vector< vector<float> > & modesshape, const vector<float> & se, \
						const vector<float> & ishape, const vector< vector<float> > & modesint, const vector<float> & ie, \
						const int & M, const vector<float> & errs, const vector<short> & vmaskin)
{
	smean=mshape;
	smodes=modesshape;
	seigs=se;
	imean=ishape;
	imodes=modesint;
	ieigs=ie;
	NumberOfSubjects=M;
	Errs=errs;
	stmask=vmaskin;
	MODE_FOUND=false;
	mode=0;
	
}



std::vector<float> shapeModel::getDeformedGrid( const std::vector<float>  & vars ) const 
{
	std::vector<float> newshape=smean;
	std::vector<float>::const_iterator sqrtseigs_i=sqrtseigs.begin();
	std::vector< std::vector<float> >::const_iterator smodes_i = smodes.begin();
	for (std::vector<float>::const_iterator vars_i=vars.begin(); vars_i!=vars.end(); vars_i++,sqrtseigs_i++, smodes_i++)
	{
		std::vector<float>::iterator new_i=newshape.begin();
		for (std::vector<float>::const_iterator smodes_j= smodes_i->begin(); smodes_j!= smodes_i->end(); smodes_j++,new_i++)
			*new_i += (*vars_i) * (*sqrtseigs_i) * (*smodes_j);//* 
	}
	
	return newshape;
}

vector<float> shapeModel::getDeformedIGrid( const vector<float> & vars) const {
	vector<float> newigrid=imean;
	for (unsigned int i=0; i< vars.size();i++){
		for (unsigned int j=0; j< imean.size(); j++){
			newigrid.at(j)+=vars.at(i)*sqrtseigs.at(i)*imodes.at(i).at(j);
		}
	}
	return newigrid;
}

	void shapeModel::registerModel(const vector< vector<float> > & flirtmat)
	{
		//start by registering mesh 
		float f_11= flirtmat.at(0).at(0);
		float f_12= flirtmat.at(0).at(1);
		float f_13= flirtmat.at(0).at(2);
		float f_14= flirtmat.at(0).at(3);
		float f_21= flirtmat.at(1).at(0);
		float f_22= flirtmat.at(1).at(1);
		float f_23= flirtmat.at(1).at(2);
		float f_24= flirtmat.at(1).at(3);
		float f_31= flirtmat.at(2).at(0);
		float f_32= flirtmat.at(2).at(1);
		float f_33= flirtmat.at(2).at(2);
		float f_34= flirtmat.at(2).at(3);
		//assumes last line is 0 0 0 1
		for (vector<float>::iterator i=smean.begin();i!=smean.end();i+=3)
		{
			float x= f_11 * (*i)  + f_12 * (*(i+1)) + f_13 * (*(i+2)) +  f_14;
			float y= f_21 * (*i)  + f_22 * (*(i+1)) + f_23 * (*(i+2)) +  f_24;
			float z= f_31 * (*i)  + f_32 * (*(i+1)) + f_33 * (*(i+2)) +  f_34;
		
			(*i)=x;
			(*(i+1))=y;
			(*(i+2))=z;				
		}
		
		for (unsigned int i=0; i<smodes.size();i++){
			for (unsigned int j=0; j<smodes.at(0).size();j+=3){
			//exclude translation from mode of variation so that it is compatible with implementation	
				float x= flirtmat.at(0).at(0)*smodes.at(i).at(j)  + flirtmat.at(0).at(1)*smodes.at(i).at(j+1) + flirtmat.at(0).at(2)*smodes.at(i).at(j+2);// +  flirtmat.at(0).at(3);
				float y= flirtmat.at(1).at(0)*smodes.at(i).at(j)  + flirtmat.at(1).at(1)*smodes.at(i).at(j+1) + flirtmat.at(1).at(2)*smodes.at(i).at(j+2) ;//+  flirtmat.at(1).at(3);
				float z= flirtmat.at(2).at(0)*smodes.at(i).at(j)  + flirtmat.at(2).at(1)*smodes.at(i).at(j+1) + flirtmat.at(2).at(2)*smodes.at(i).at(j+2) ;//+  flirtmat.at(2).at(3);
		
				smodes.at(i).at(j)=x;
				smodes.at(i).at(j+1)=y;
				smodes.at(i).at(j+2)=z;	
			}
			
		}

	}
}
