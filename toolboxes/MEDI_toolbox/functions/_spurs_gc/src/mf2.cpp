//--------------------------------------------------------------
// file: mf2.cpp - Matlab Mex function application
//
// Debugging notes:
//   (1) In the Project Settings, under the Debug tab, set
//       the "Executable for debug session" to:
//         C:\MATLAB6p1\bin\win32\matlab.exe
//       and set the working directory to:
//         debug
//       and (optionally) set the "Program arguments" to:
//          /nosplash
//   (2) From the menu bar, select Debug->Exceptions to bring up 
//       the exception handling dialog. Change the trapping of 
//       "access violation" from "Stop in not handled" to "Stop
//       always." This will allow the lines of code causing
//       exceptions to be readily identified.
//   (3) When you begin a debug session, Matlab is invoked.
//       Your mex code will be executed when the mex function
//       is called from Matlab.
//   (4) You may use a simple startup.m file to call your mex 
//       function as soon as Matlab starts
//   (5) A virus scanner may cause long delays when Matlab is
//       starting. It may be worthwhile to temporarily disable
//       the virus scanner while you are developing.
// Development Notes:
//   (1) All global variables in your mex function remain in
//       scope after the mex function is executed the first time
//       in a Matlab session.
//   (2) Likewise, all static variables are persistent between
//       function calls.
//   (3) The "clear functions" command may be used to unload
//       your mex function DLL, and force all globals to be
//       reinitialized the next time the mex function is called.
//   (4) Do not use malloc() or free(). Use mxCalloc() and 
//       mxFree() instead. 
//   (5) If you are writing code to be portable to other
//       environments, you may use the MATLAB_MEX_FILE macro
//       to determine if the code is targeting a mex file.
//--------------------------------------------------------------
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include <string.h>
#include "print.h"

extern "C" {
#include "mex.h"
}

// TODO: Add you supporting functions here


//--------------------------------------------------------------
// function: mf2 - Entry point from Matlab environment (via 
//   mexFucntion(), below)
// INPUTS:
//   nlhs - number of left hand side arguments (outputs)
//   plhs[] - pointer to table where created matrix pointers are
//            to be placed
//   nrhs - number of right hand side arguments (inputs)
//   prhs[] - pointer to table of input matrices
//--------------------------------------------------------------
void mf2( int nlhs, mxArray *plhs[], int nrhs, const mxArray  *prhs[] )
{
  
  //Declarations
  float *sourcesinkVal;
  float *remainVal;
  float *outArrayFlow,*outArrayAssign;
  mwSize dims[1] = {1};
  
  mwSize sourcerowLen, sourcecolLen;
  mwSize remainrowLen, remaincolLen;
  
  //Get matrix sourcesink
  sourcesinkVal = (float *) mxGetPr(prhs[0]);
  sourcerowLen = (mwSize) mxGetN(prhs[0]);
  sourcecolLen = (mwSize) mxGetM(prhs[0]);
  
  //for(int j=0; j<10; j++) {
  //    for (int i=0; i<3; i++) {
  //        PRINT(sourcesinkVal[i*sourcecolLen+j]);
  //    }
  //}
  
  //PRINT(sourcerowLen);
  //PRINT(sourcecolLen);
	
  //Get matrix remain
  remainVal    = (float *) mxGetPr(prhs[1]);
  remainrowLen = (mwSize) mxGetN(prhs[1]);
  remaincolLen = (mwSize) mxGetM(prhs[1]);	

  //PRINT(remainrowLen);
  //PRINT(remaincolLen);
  
  Graph::node_id *nodes = new Graph::node_id[sourcecolLen];	// dynamic allocation
  //PRINT(__LINE__);

  Graph *g = new Graph();
  //PRINT(__LINE__);
	
  for(int t=0; t<sourcecolLen; t++)
      //nodes[t] = g -> add_node();
	  nodes[(int)sourcesinkVal[t]-1] = g -> add_node();
  
  //PRINT(__LINE__);

  /* capacity of arcs between terminals and nodes*/
  for(int k=0; k<(sourcecolLen); k++) // sourcesink
      g -> set_tweights(nodes[(int)sourcesinkVal[k]-1], sourcesinkVal[sourcecolLen + k], sourcesinkVal[2*sourcecolLen + k]); 
  //PRINT(__LINE__);
  	
  int  inc=0;
  for(int h=0; h<remaincolLen; h++)	
		g -> add_edge(nodes[(int)remainVal[h]-1],nodes[(int)remainVal[remaincolLen + h]-1], remainVal[2*remaincolLen + h],remainVal[3*remaincolLen + h]);
  //PRINT(__LINE__);
  
  Graph::flowtype flow = g -> maxflow();

  //PRINT(__LINE__);
  
  //Allocate memory and assign output pointer
  //plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type
  plhs[0] = mxCreateNumericArray(1,dims, mxSINGLE_CLASS, mxREAL);
  mwSize dims2[] = {sourcecolLen,2};
  plhs[1] = mxCreateNumericArray(2,dims2, mxSINGLE_CLASS, mxREAL);
  
  //PRINT(__LINE__);
  
  //Get a pointer to the data space in our newly allocated memory
  outArrayFlow = (float*) mxGetPr(plhs[0]);
  //PRINT(__LINE__);
  
  outArrayFlow[0] = flow;
  //PRINT(__LINE__);
  
  //Get a pointer to the data space in our newly allocated memory
  outArrayAssign = (float*) mxGetPr(plhs[1]);
  //PRINT(__LINE__);

  for(int i=0;i<sourcecolLen;i++)
{
    
	//assigns the number of the nodes
	outArrayAssign[i] = sourcesinkVal[i];
	
    //If o n� becomes assigned to source we assign the value 0
	//Otherwise we assign the value 1.
	if (g->what_segment(nodes[(int)sourcesinkVal[i]-1]) == Graph::SOURCE)
		
		outArrayAssign[sourcecolLen + i] = 0;
			
	else
		outArrayAssign[sourcecolLen + i] = 1;
	
}
    //PRINT(__LINE__);

  delete g;
    //PRINT(__LINE__);


} // end mf2()

extern "C" {
  //--------------------------------------------------------------
  // mexFunction - Entry point from Matlab. From this C function,
  //   simply call the C++ application function, above.
  //--------------------------------------------------------------
  void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray  *prhs[] )
  {
    mf2(nlhs, plhs, nrhs, prhs);
  }



}


