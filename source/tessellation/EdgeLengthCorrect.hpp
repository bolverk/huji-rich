#include "tessellation.hpp"

/*
\brief Corrects the length of the edges to be second order in time
\param tessold The tessellation at the start of the time step
\param tessnew The tessellation at the end of the time step
\param lengths This is the output of the correct edges lengths
*/
void CorrectEdgeLength(Tessellation const& tessold,Tessellation const& tessnew,
	vector<double> &lengths);

/*
\brief Corrects the length of the edges to be second order in time
\param tessold The tessellation at the start of the time step
\param tessold The tessellation at the middle of the time step
\param tessnew The tessellation at the end of the time step
\param lengths This is the output of the correct edges lengths
*/
void CorrectEdgeLength(Tessellation const& tessold,Tessellation const& tessmid,
	Tessellation const& tessnew,vector<double> &lengths);