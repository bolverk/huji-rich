#include "mpi_macro.hpp"
#ifdef RICH_MPI

int get_mpi_rank(void)
{
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int get_mpi_size(void)
{
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

int MPI_Send_Vector2D(Vector2D const& vec,int dest,int tag, MPI_Comm comm)
{
  double temp[2];
  temp[0]=vec.x;
  temp[1]=vec.y;
  int err=MPI_Send(&temp[0],2,MPI_DOUBLE,dest,tag,comm);
  return err;
}

int MPI_VectorSend_Vector2D(vector<Vector2D> const& vec,int dest,int tag, MPI_Comm comm)
{
  int n=(int)vec.size();
  if(vec.empty())
    {
      double temp=0;
      int err=MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
      return err;
    }
  vector<double> temp(n*2);
  for(int i=0;i<n;++i)
    {
      temp[2*i]=vec[i].x;
      temp[2*i+1]=vec[i].y;
    }
  int err;
  err=MPI_Send(&temp[0],2*n,MPI_DOUBLE,dest,tag,comm);
  return err;
}

int MPI_Recv_Vector2D(Vector2D &vec,int source, int tag, MPI_Comm comm)
{
  double temp[2];
  int err=0;
  err=MPI_Recv(&temp[0],2,MPI_DOUBLE,source,tag,comm,MPI_STATUS_IGNORE);
  vec.x=temp[0];
  vec.y=temp[1];
  return err;
}

int MPI_VectorRecv_Vector2D(vector<Vector2D> &vec,int source, int tag, MPI_Comm comm)
{
  vec.clear();
  MPI_Status status;
  MPI_Probe(source,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
  if(status.MPI_TAG==1)
    {
      double temp;
      int err=MPI_Recv(&temp,1,MPI_DOUBLE,source,1,comm,&status);
      return err;
    }
  int n=0;
  MPI_Get_count(&status,MPI_DOUBLE,&n);
  vector<double> temp(n);
  vec.resize(n/2);
  n/=2;
  int err=0;
  err=MPI_Recv(&temp[0],2*n,MPI_DOUBLE,source,tag,comm,MPI_STATUS_IGNORE);
  for(int i=0;i<n;++i)
    {
      vec[i].x=temp[2*i];
      vec[i].y=temp[2*i+1];
    }
  return err;
}

int MPI_VectorBcast_Vector2D(vector<Vector2D> &vec,int root, MPI_Comm comm,int rank)
{
  int n=(int)vec.size();
  vector<double> temp(n*2);
  if(rank==root)
    {
      for(int i=0;i<n;++i)
	{
	  temp[2*i]=vec[i].x;
	  temp[2*i+1]=vec[i].y;
	}
    }
  int err=MPI_Bcast(&temp[0],n*2,MPI_DOUBLE,root,comm);
  if(rank!=root)
    {
      for(int i=0;i<n;++i)
	{
	  vec[i].x=temp[2*i];
	  vec[i].y=temp[2*i+1];
	}
    }
  return err;
}

bool PointInsideCell(Tessellation const& tess,int cell_index,Vector2D const & point)
{
  boost::array<Vector2D,3> triangle;
  vector<int>const& cell_edges=tess.GetCellEdges(cell_index);
  triangle[0]=tess.GetEdge(cell_edges[0]).GetVertex(0);
  triangle[1]=tess.GetEdge(cell_edges[0]).GetVertex(1);
  triangle[2]=tess.GetCellCM(cell_index);
  double temp=orient2d(triangle);
  triangle[2]=point;
  return (orient2d(triangle)*temp)>0;
}

/*void BuildTree(ANNkd_tree *&tree,ANNpointArray &treePoints,
  Tessellation const& tess)
  {
  int N=tess.GetPointNo();
  treePoints=annAllocPts(N,2);
  for(int i=0;i<N;++i)
  {
  treePoints[i][0]=tess.GetMeshPoint(i).x;
  treePoints[i][1]=tess.GetMeshPoint(i).y;
  }
  tree=new ANNkd_tree(treePoints,N,2,1,ANN_KD_SUGGEST);
  return;
  }

  void KillTree(ANNkd_tree *&tree,ANNpointArray &treePoints)
  {
  annDeallocPts(treePoints);
  delete tree;
  annClose();	
  }

  int FindContainingCell(ANNkd_tree *tree,Vector2D const& point)
  {
  ANNpoint queryPt;
  queryPt=annAllocPt(2);
  queryPt[0]=point.x;
  queryPt[1]=point.y;
  ANNidxArray nnIdx; // near neighbor indices
  ANNdistArray dists; // near neighbor distances
  nnIdx = new ANNidx[1]; // allocate near neigh indices
  dists = new ANNdist[1]; // allocate near neighbor dists
  tree->annkSearch(queryPt,1,nnIdx,dists);
  int res=nnIdx[0];
  delete [] nnIdx;
  delete [] dists;
  annDeallocPt(queryPt);
  return res;
  }
*/
void ConvertVector2DToDouble(vector<Vector2D> const& vec,vector<double> &res)
{
  int n=(int) vec.size();
  res.resize(n*2);
  for(int i=0;i<n;++i)
    {
      res[2*i]=vec[i].x;
      res[2*i+1]=vec[i].y;
    }
}

vector<Vector2D> MPI_MassSendRecvVectorVector2D
(vector<vector<Vector2D> > const& tosend,
 vector<int> const& proclist,vector<int> const& procorder)
{
  vector<Vector2D> res;
  int n=(int)procorder.size();
  int nlist=(int)proclist.size();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  for(int i=0;i<n;++i)
    {
      int index=find(proclist.begin(),proclist.end(),procorder[i])-proclist.begin();
      // Do we talk with this processor?
      if(index<nlist)
	{
	  if(rank<procorder[i])
	    {
	      if(tosend[index].empty())
		{
		  double temp=0;
		  MPI_Send(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD);
		}
	      else
		MPI_VectorSend_Vector2D(tosend[index],procorder[i],0,
					MPI_COMM_WORLD);
	      MPI_Status status;
	      MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	      if(status.MPI_TAG==0)
		{
		  vector<Vector2D> tempres;
		  MPI_VectorRecv_Vector2D(tempres,procorder[i],0,MPI_COMM_WORLD);
		  res.insert(res.end(),tempres.begin(),tempres.end());
		}
	      else
		{
		  double temp;
		  MPI_Recv(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,&status);
		}
	    }
	  else
	    {
	      MPI_Status status;
	      MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	      if(status.MPI_TAG==0)
		{
		  vector<Vector2D> tempres;
		  MPI_VectorRecv_Vector2D(tempres,procorder[i],0,MPI_COMM_WORLD);
		  res.insert(res.end(),tempres.begin(),tempres.end());
		}
	      else
		{
		  double temp;
		  MPI_Recv(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,&status);
		}
	      if(tosend[index].empty())
		{
		  double temp=0;
		  MPI_Send(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD);
		}
	      else
		MPI_VectorSend_Vector2D(tosend[index],procorder[i],0,
					MPI_COMM_WORLD);
	    }
	}
    }
  return res;
}

vector<int> GetProcOrder(int rank,int worldsize)
{
  if(worldsize==1)
    return vector<int> ();
  vector<int> procorder(worldsize-1);
  if(rank==0)
    {
      for(int i=0;i<(worldsize-1);++i)
	procorder[i]=i+1;
    }
  if(rank!=(worldsize-1))
    {
      for(int i=0;i<(worldsize-1);++i)
	{
	  int temp=(i-rank+worldsize)%(worldsize-1);
	  if(temp!=rank)
	    procorder[i]=temp;
	  else
	    procorder[i]=worldsize-1;
	}
    }
  else
    {
      if(worldsize>1)
	{
	  int half=(int)(worldsize/2);
	  for(int i=0;i<(worldsize-1);++i)
	    {
	      if(i%2==0)
		procorder[i]=half+i/2;
	      else
		procorder[i]=(int)(i/2)+1;
	    }
	  procorder[worldsize-2]=0;
	}
    }
  return procorder;
}

void SendRecvExtensive(vector<Conserved> const& cons,vector<vector<double> > const&
		       tracers,vector<size_t> const& customevolutions,vector<vector<int> > const& sentcells,
		       vector<int> const& sentprocs,vector<Conserved> &ptoadd,vector<vector<double> >
		       &ttoadd,vector<size_t> &ctoadd)
{
  bool traceractive=tracers.empty() ? false : true;
  if(traceractive)
    traceractive=tracers[0].empty() ? false : true;

  // Take care of self send hydro
  int nlist=(int)sentprocs.size();
  ptoadd.clear();
  ttoadd.clear();
  ctoadd.clear();

  // Talk with other procs
  int rank,worldsize;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&worldsize);
  const vector<int> procorder=GetProcOrder(rank,worldsize);

  int n=worldsize-1;
  for(int i=0;i<n;++i)
    {
      int index=find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
      // Do we talk with this processor?
      if(index<nlist)
	{
	  if(rank<procorder[i])
	    {
	      vector<Conserved> ptemp=VectorValues(cons,sentcells[index]);
	      MPI_SendVectorConserved(ptemp,procorder[i],0,MPI_COMM_WORLD);
	      MPI_RecvVectorConserved(ptemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD);
	      if(!ptemp.empty())
		{
		  ptoadd.insert(ptoadd.end(),ptemp.begin(),ptemp.end());
		}
	      if(traceractive)
		{
		  vector<vector<double> > ttemp=VectorValues(tracers,sentcells[index]);
		  MPI_SendVectorTracer(ttemp,procorder[i],0,MPI_COMM_WORLD);
		  MPI_RecvVectorTracer(ttemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,
				       (int)tracers[0].size());
		  ttoadd.insert(ttoadd.end(),ttemp.begin(),ttemp.end());
		}
	      if(!sentcells[index].empty())
		{
		  vector<size_t> ctemp=VectorValues(customevolutions,sentcells[index]);
		  MPI_Send(&ctemp[0],(int)sentcells[index].size(),MPI_UNSIGNED,procorder[i],0,
			   MPI_COMM_WORLD);
		}
	      else
		{
		  int temp=0;
		  MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
		}
	      MPI_Status status;
	      MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	      if(status.MPI_TAG==1)
		{
		  int temp;
		  MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,&status);
		}
	      else
		{
		  int count;
		  MPI_Get_count(&status,MPI_UNSIGNED,&count);
		  vector<size_t> ctemp(count);
		  MPI_Recv(&ctemp[0],count,MPI_UNSIGNED,procorder[i],status.MPI_TAG,
			   MPI_COMM_WORLD,&status);
		  ctoadd.insert(ctoadd.end(),ctemp.begin(),ctemp.end());
		}
	    }
	  else
	    {
	      vector<Conserved> ptemp;
	      MPI_RecvVectorConserved(ptemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD);
	      if(!ptemp.empty())
		{
		  ptoadd.insert(ptoadd.end(),ptemp.begin(),ptemp.end());
		}
	      ptemp=VectorValues(cons,sentcells[index]);
	      MPI_SendVectorConserved(ptemp,procorder[i],0,MPI_COMM_WORLD);
	      if(traceractive)
		{
		  vector<vector<double> > ttemp;
		  MPI_RecvVectorTracer(ttemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,
				       (int)tracers[0].size());
		  ttoadd.insert(ttoadd.end(),ttemp.begin(),ttemp.end());
		  ttemp=VectorValues(tracers,sentcells[index]);
		  MPI_SendVectorTracer(ttemp,procorder[i],0,MPI_COMM_WORLD);
		}
	      MPI_Status status;
	      MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	      if(status.MPI_TAG==1)
		{
		  int temp;
		  MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,&status);
		}
	      else
		{
		  int count;
		  MPI_Get_count(&status,MPI_UNSIGNED,&count);
		  vector<size_t> ctemp(count);
		  MPI_Recv(&ctemp[0],count,MPI_UNSIGNED,procorder[i],status.MPI_TAG,
			   MPI_COMM_WORLD,&status);
		  ctoadd.insert(ctoadd.end(),ctemp.begin(),ctemp.end());
		}
	      if(!sentcells[index].empty())
		{
		  vector<size_t> ctemp=VectorValues(customevolutions,sentcells[index]);
		  MPI_Send(&ctemp[0],(int)sentcells[index].size(),MPI_UNSIGNED,procorder[i],0,
			   MPI_COMM_WORLD);
		}
	      else
		{
		  int temp=0;
		  MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
		}
	    }
	}
    }
}

void SendRecvShockedStatus(vector<char> const& shockedcells,
			   vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
			   vector<char> &btoadd)
{
  // Take care of self send cells
  int nlist=(int)sentprocs.size();
  btoadd.clear();
	
  // Talk with other procs
  int rank,worldsize;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&worldsize);
  const vector<int> procorder=GetProcOrder(rank,worldsize);

  int n=worldsize-1;
  for(int i=0;i<n;++i)
    {
      int index=find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
      // Do we talk with this processor?
      if(index<nlist)
	{
	  if(rank<procorder[i])
	    {
	      vector<char> btemp=VectorValues(shockedcells,sentcells[index]);
	      if(!btemp.empty())
		MPI_Send(&btemp[0],(int)btemp.size(),MPI_CHAR,procorder[i],0,
			 MPI_COMM_WORLD);
	      else
		{
		  char bbtemp=0;
		  MPI_Send(&bbtemp,1,MPI_CHAR,procorder[i],1,
			   MPI_COMM_WORLD);
		}
	      MPI_Status status;
	      MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	      if(status.MPI_TAG==1)
		{
		  char bbtemp;
		  MPI_Recv(&bbtemp,1,MPI_CHAR,procorder[i],1,MPI_COMM_WORLD,
			   &status);
		}
	      else
		{
		  int count;
		  MPI_Get_count(&status,MPI_CHAR,&count);
		  vector<char> toadd(count);
		  MPI_Recv(&toadd[0],count,MPI_CHAR,procorder[i],0,
			   MPI_COMM_WORLD,&status);
		  btoadd.insert(btoadd.end(),toadd.begin(),toadd.end());
		}
	    }
	  else
	    {
	      MPI_Status status;
	      MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	      if(status.MPI_TAG==1)
		{
		  char bbtemp;
		  MPI_Recv(&bbtemp,1,MPI_CHAR,procorder[i],1,MPI_COMM_WORLD,
			   &status);
		}
	      else
		{
		  int count;
		  MPI_Get_count(&status,MPI_CHAR,&count);
		  vector<char> toadd(count);
		  MPI_Recv(&toadd[0],count,MPI_CHAR,procorder[i],0,
			   MPI_COMM_WORLD,&status);
		  btoadd.insert(btoadd.end(),toadd.begin(),toadd.end());
		}
	      vector<char> btemp=VectorValues(shockedcells,sentcells[index]);
	      if(!btemp.empty())
		MPI_Send(&btemp[0],(int)btemp.size(),MPI_CHAR,procorder[i],0,
			 MPI_COMM_WORLD);
	      else
		{
		  char bbtemp=0;
		  MPI_Send(&bbtemp,1,MPI_CHAR,procorder[i],1,
			   MPI_COMM_WORLD);
		}
	    }
	}
    }
}

void SendRecvVectorDouble(vector<double> const& vec,
			  vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
			  vector<double> &toadd)
{
  // Take care of self send cells
  int nlist=(int)sentprocs.size();
  toadd.clear();
	
  // Talk with other procs
  int rank,worldsize;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&worldsize);
  const vector<int> procorder=GetProcOrder(rank,worldsize);

  int n=worldsize-1;
  for(int i=0;i<n;++i)
    {
      int index=find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
      // Do we talk with this processor?
      if(index<nlist)
	{
	  if(rank<procorder[i])
	    {
	      vector<double> temp=VectorValues(vec,sentcells[index]);
	      if(!temp.empty())
		MPI_Send(&temp[0],(int)temp.size(),MPI_DOUBLE,procorder[i],0,
			 MPI_COMM_WORLD);
	      else
		{
		  double bbtemp=0;
		  MPI_Send(&bbtemp,1,MPI_DOUBLE,procorder[i],1,
			   MPI_COMM_WORLD);
		}
	      MPI_Status status;
	      MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	      if(status.MPI_TAG==1)
		{
		  double bbtemp;
		  MPI_Recv(&bbtemp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,
			   &status);
		}
	      else
		{
		  int count;
		  MPI_Get_count(&status,MPI_DOUBLE,&count);
		  vector<double> ttoadd(count);
		  MPI_Recv(&ttoadd[0],count,MPI_DOUBLE,procorder[i],0,
			   MPI_COMM_WORLD,&status);
		  toadd.insert(toadd.end(),ttoadd.begin(),ttoadd.end());
		}
	    }
	  else
	    {
	      MPI_Status status;
	      MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	      if(status.MPI_TAG==1)
		{
		  double bbtemp;
		  MPI_Recv(&bbtemp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,
			   &status);
		}
	      else
		{
		  int count;
		  MPI_Get_count(&status,MPI_DOUBLE,&count);
		  vector<double> ttoadd(count);
		  MPI_Recv(&ttoadd[0],count,MPI_DOUBLE,procorder[i],0,
			   MPI_COMM_WORLD,&status);
		  toadd.insert(toadd.end(),ttoadd.begin(),ttoadd.end());
		}
	      vector<double> temp=VectorValues(vec,sentcells[index]);
	      if(!temp.empty())
		MPI_Send(&temp[0],(int)temp.size(),MPI_DOUBLE,procorder[i],0,
			 MPI_COMM_WORLD);
	      else
		{
		  double bbtemp=0;
		  MPI_Send(&bbtemp,1,MPI_DOUBLE,procorder[i],1,
			   MPI_COMM_WORLD);
		}
	    }
	}
    }
}


void SendRecvOldVector2D(vector<Vector2D> const& points,
			 vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
			 vector<Vector2D> &toadd)
{

  int nlist=(int)sentprocs.size();
  toadd.clear();

  // Talk with other procs
  int rank,worldsize;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&worldsize);
  const vector<int> procorder=GetProcOrder(rank,worldsize);

  int n=worldsize-1;
  for(int i=0;i<n;++i)
    {
      int index=find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
      // Do we talk with this processor?
      if(index<nlist)
	{
	  if(rank<procorder[i])
	    {
	      vector<Vector2D> vtemp=VectorValues(points,sentcells[index]);
	      MPI_VectorSend_Vector2D(vtemp,procorder[i],0,MPI_COMM_WORLD);
	      MPI_VectorRecv_Vector2D(vtemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD);
	      if(!vtemp.empty())
		{
		  toadd.insert(toadd.end(),vtemp.begin(),vtemp.end());
		}
	    }
	  else
	    {
	      vector<Vector2D> vtemp;
	      MPI_VectorRecv_Vector2D(vtemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD);
	      if(!vtemp.empty())
		{
		  toadd.insert(toadd.end(),vtemp.begin(),vtemp.end());
		}
	      vtemp=VectorValues(points,sentcells[index]);
	      MPI_VectorSend_Vector2D(vtemp,procorder[i],0,MPI_COMM_WORLD);
	    }				
	}
    }
}


void KeepLocalPoints(vector<Conserved> &cons,vector<vector<double> > &tracers,
		     vector<size_t> &customevolutions,vector<int> const& localpoints)
{
  bool traceractive=tracers.empty() ? false : true;
  if(traceractive)
    traceractive=tracers[0].empty() ? false : true;
  cons=VectorValues(cons,localpoints);
  if(traceractive)
    tracers=VectorValues(tracers,localpoints);
  customevolutions=VectorValues(customevolutions,localpoints);
}

void SendRecvHydro(vector<Primitive> &cells,vector<vector<double> > &tracers,
		   vector<size_t> &customevolutions,vector<vector<int> > sentcells,
		   vector<int> sentprocs,EquationOfState const& eos,int totalpoints)
{
  bool traceractive=tracers.empty() ? false : true;
  if(traceractive)
    traceractive=tracers[0].empty() ? false : true;
  // Remove all the sent points
  /*vector<Primitive> rescells=VectorValues(cells,selfindex);
    vector<vector<double> > restracer=VectorValues(tracers,selfindex);
    vector<size_t> rescustomevolutions=VectorValues(customevolutions,selfindex);
  */
  // Take care of self send hydro
  int nlist=(int)sentprocs.size();
  vector<Primitive> ptoadd;
  vector<vector<double> > ttoadd;
  vector<size_t> ctoadd;
  for(int i=0;i<nlist;++i)
    {
      if(sentprocs[i]==-1)
	{
	  if(sentcells[i].empty())
	    continue;
	  vector<Primitive> ptemp=VectorValues(cells,sentcells[i]);
	  ptoadd.insert(ptoadd.end(),ptemp.begin(),ptemp.end());
	  vector<vector<double> > ttemp=VectorValues(tracers,sentcells[i]);
	  if(traceractive)
	    ttoadd.insert(ttoadd.end(),ttemp.begin(),ttemp.end());
	  vector<size_t> ctemp=VectorValues(customevolutions,sentcells[i]);
	  ctoadd.insert(ctoadd.end(),ctemp.begin(),ctemp.end());
	}
    }

  // Talk with other procs
  int rank,worldsize;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&worldsize);
  const vector<int> procorder=GetProcOrder(rank,worldsize);
  int n=worldsize-1;
  vector<int> sentprocs2(sentprocs);
  vector<int> talkedto;
  vector<vector<Primitive> > pmany(sentprocs.size());
  vector<vector<vector<double> > > tmany(sentprocs.size());
  vector<vector<size_t> > cmany(sentprocs.size());
  for(int k=0;k<3;++k)
    {
      for(int i=0;i<n;++i)
	{
	  nlist=(int)sentprocs.size();
	  int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
	  // Do we talk with this processor?
	  if(index<nlist)
	    {
	      if(rank<procorder[i])
		{
		  vector<Primitive> ptemp=VectorValues(cells,sentcells[index]);
		  MPI_SendVectorPrimitive(ptemp,procorder[i],0,MPI_COMM_WORLD);
		  MPI_RecvVectorPrimitive(ptemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,
					  eos);
		  if(!ptemp.empty())
		    {
		      //ptoadd.insert(ptoadd.end(),ptemp.begin(),ptemp.end());
		      pmany[talkedto.size()]=ptemp;
		    }
		  if(traceractive)
		    {
		      vector<vector<double> > ttemp=VectorValues(tracers,sentcells[index]);
		      MPI_SendVectorTracer(ttemp,procorder[i],0,MPI_COMM_WORLD);
		      MPI_RecvVectorTracer(ttemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,
					   (int)tracers[0].size());
		      //ttoadd.insert(ttoadd.end(),ttemp.begin(),ttemp.end());
		      tmany[talkedto.size()]=ttemp;
		    }
		  if(!sentcells[index].empty())
		    {
		      vector<size_t> ctemp=VectorValues(customevolutions,sentcells[index]);
		      MPI_Send(&ctemp[0],(int)sentcells[index].size(),MPI_UNSIGNED,procorder[i],0,
			       MPI_COMM_WORLD);
		    }
		  else
		    {
		      int temp=0;
		      MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
		    }
		  MPI_Status status;
		  MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		  if(status.MPI_TAG==1)
		    {
		      int temp;
		      MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,&status);
		    }
		  else
		    {
		      int count;
		      MPI_Get_count(&status,MPI_UNSIGNED,&count);
		      vector<size_t> ctemp(count);
		      MPI_Recv(&ctemp[0],count,MPI_UNSIGNED,procorder[i],status.MPI_TAG,
			       MPI_COMM_WORLD,&status);
		      //ctoadd.insert(ctoadd.end(),ctemp.begin(),ctemp.end());
		      cmany[talkedto.size()]=ctemp;
		    }
		}
	      else
		{
		  vector<Primitive> ptemp;
		  MPI_RecvVectorPrimitive(ptemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,
					  eos);
		  if(!ptemp.empty())
		    {
		      //ptoadd.insert(ptoadd.end(),ptemp.begin(),ptemp.end());
		      pmany[talkedto.size()]=ptemp;
		    }
		  ptemp=VectorValues(cells,sentcells[index]);
		  MPI_SendVectorPrimitive(ptemp,procorder[i],0,MPI_COMM_WORLD);
		  if(traceractive)
		    {
		      vector<vector<double> > ttemp;
		      MPI_RecvVectorTracer(ttemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,
					   (int)tracers[0].size());
		      //ttoadd.insert(ttoadd.end(),ttemp.begin(),ttemp.end());
		      tmany[talkedto.size()]=ttemp;
		      ttemp=VectorValues(tracers,sentcells[index]);
		      MPI_SendVectorTracer(ttemp,procorder[i],0,MPI_COMM_WORLD);
		    }
		  MPI_Status status;
		  MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		  if(status.MPI_TAG==1)
		    {
		      int temp;
		      MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,&status);
		    }
		  else
		    {
		      int count;
		      MPI_Get_count(&status,MPI_UNSIGNED,&count);
		      vector<size_t> ctemp(count);
		      MPI_Recv(&ctemp[0],count,MPI_UNSIGNED,procorder[i],status.MPI_TAG,
			       MPI_COMM_WORLD,&status);
		      //ctoadd.insert(ctoadd.end(),ctemp.begin(),ctemp.end());
		      cmany[talkedto.size()]=ctemp;
		    }
		  if(!sentcells[index].empty())
		    {
		      vector<size_t> ctemp=VectorValues(customevolutions,sentcells[index]);
		      MPI_Send(&ctemp[0],(int)sentcells[index].size(),MPI_UNSIGNED,procorder[i],0,
			       MPI_COMM_WORLD);
		    }
		  else
		    {
		      int temp=0;
		      MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
		    }
		}
	      sentprocs.erase(sentprocs.begin()+index);
	      sentcells.erase(sentcells.begin()+index);
	      talkedto.push_back(procorder[i]);
	    }
	}
    }
  if(!pmany.empty())
    {
      int nn=(int)sentprocs2.size();
      for(int i=0;i<nn;++i)
	{
	  //	int index=Find(talkedto.begin(),talkedto.end(),sentprocs2[i])-talkedto.begin();
	  if(!pmany[i].empty())
	    {
	      ptoadd.insert(ptoadd.end(),pmany[i].begin(),pmany[i].end());
	      ctoadd.insert(ctoadd.end(),cmany[i].begin(),cmany[i].end());
	      if(traceractive)
		ttoadd.insert(ttoadd.end(),tmany[i].begin(),tmany[i].end());
	    }
	  //talkedto.erase(talkedto.begin()+index);
	  //pmany.erase(pmany.begin()+index);
	  //cmany.erase(cmany.begin()+index);
	  //tmany.erase(tmany.begin()+index);
	}
    }
  if(!ptoadd.empty())
    {
      cells.resize((size_t)totalpoints-ptoadd.size());
      cells.insert(cells.end(),ptoadd.begin(),ptoadd.end());
    }
  //cells=rescells;
  if(!ttoadd.empty())
    {
      tracers.resize((size_t)totalpoints-ttoadd.size());
      tracers.insert(tracers.end(),ttoadd.begin(),ttoadd.end());
    }
  //tracers=restracer;
  if(!ctoadd.empty())
    {
      customevolutions.resize((size_t)totalpoints-ctoadd.size());
      customevolutions.insert(customevolutions.end(),ctoadd.begin(),
			      ctoadd.end());
    }
  //customevolutions=rescustomevolutions;
}

int MPI_SendVectorPrimitive(vector<Primitive> const& vec,int dest,int tag,
			    MPI_Comm comm)
{
  if(vec.empty())
    {
      double temp=0;
      return MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
    }
  int n=(int)vec.size();
  vector<double> tosend(n*4);
  for(int i=0;i<n;++i)
    {
      tosend[4*i]=vec[i].Density;
      tosend[4*i+1]=vec[i].Pressure;
      tosend[4*i+2]=vec[i].Velocity.x;
      tosend[4*i+3]=vec[i].Velocity.y;
    }
  return MPI_Send(&tosend[0],4*n,MPI_DOUBLE,dest,tag,comm);
}

void MPI_RecvVectorPrimitive(vector<Primitive> &vec,int dest,int tag,
			     MPI_Comm comm,EquationOfState const& eos)
{
  vec.clear();
  MPI_Status status;
  MPI_Probe(dest,tag,comm,&status);
  if(tag==1)
    {
      double temp;
      MPI_Recv(&temp,1,MPI_DOUBLE,dest,tag,comm,&status);
      return;
    }
  int nrecv;
  MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
  vector<double> temp(nrecv);
  MPI_Recv(&temp[0],nrecv,MPI_DOUBLE,dest,tag,comm,&status);
  int ntotal=nrecv/4;
  vec.reserve(ntotal);
  for(int i=0;i<ntotal;++i)
    {
      Primitive ptemp(temp[4*i],temp[4*i+1],Vector2D(temp[4*i+2],temp[4*i+3]),0,0);
      ptemp.Energy=eos.dp2e(ptemp.Density,ptemp.Pressure);
      ptemp.SoundSpeed=eos.dp2c(ptemp.Density,ptemp.Pressure);
      vec.push_back(ptemp);
    }
}

int MPI_SendVectorTracer(vector<vector<double> > const& vec,int dest,int tag,
			 MPI_Comm comm)
{
  if(vec.empty())
    {
      double temp=0;
      return MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
    }
  int n=(int)vec.size();
  int ntracer=(int)vec[0].size();
  vector<double> tosend(n*ntracer);
  for(int i=0;i<n;++i)
    {
      for(int j=0;j<ntracer;++j)
	tosend[ntracer*i+j]=vec[i][j];
    }
  return MPI_Send(&tosend[0],ntracer*n,MPI_DOUBLE,dest,tag,comm);
}

void MPI_RecvVectorTracer(vector<vector<double> > &vec,int dest,int tag,
			  MPI_Comm comm,int ntracer)
{
  vec.clear();
  int nrecv;
  MPI_Status status;
  MPI_Probe(dest,tag,comm,&status);
  if(status.MPI_TAG==1)
    {
      double temp;
      MPI_Recv(&temp,1,MPI_DOUBLE,dest,1,comm,&status);
      return;
    }
  MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
  vector<double> temp(nrecv);
  MPI_Recv(&temp[0],nrecv,MPI_DOUBLE,dest,MPI_ANY_TAG,comm,&status);
  int length=nrecv/ntracer;
  vec.resize(length);
  for(int i=0;i<length;++i)
    {
      vec[i].resize(ntracer);
      for(int j=0;j<ntracer;++j)
	vec[i][j]=temp[ntracer*i+j];
    }
}

int MPI_SendVectorGrad(vector<ReducedPrimitiveGradient2D> const&vec,int dest,int
		       tag,MPI_Comm comm)
{
  if(vec.empty())
    {
      double temp=0;
      return MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
    }
  int n=(int)vec.size();
  int gradlength=2*(int)vec[0].tracers.size()+8;
  vector<double> tosend(n*gradlength);
  for(int i=0;i<n;++i)
    {
      tosend[gradlength*i]=vec[i].density.x;
      tosend[gradlength*i+1]=vec[i].density.y;
      tosend[gradlength*i+2]=vec[i].pressure.x;
      tosend[gradlength*i+3]=vec[i].pressure.y;
      tosend[gradlength*i+4]=vec[i].xvelocity.x;
      tosend[gradlength*i+5]=vec[i].xvelocity.y;
      tosend[gradlength*i+6]=vec[i].yvelocity.x;
      tosend[gradlength*i+7]=vec[i].yvelocity.y;
      for(int j=0;j<(gradlength-8)/2;++j)
	{
	  tosend[gradlength*i+j+8]=vec[i].tracers[j].x;
	  tosend[gradlength*i+j+9]=vec[i].tracers[j].y;
	}
    }
  return MPI_Send(&tosend[0],gradlength*n,MPI_DOUBLE,dest,tag,comm);
}

void MPI_RecvVectorGrad(vector<ReducedPrimitiveGradient2D> &vec,int dest,int
			tag,MPI_Comm comm,int gradlength)
{
  MPI_Status status;
  MPI_Probe(dest,tag,comm,&status);
  if(status.MPI_TAG==1)
    {
      double temp;
      MPI_Recv(&temp,1,MPI_DOUBLE,dest,1,comm,&status);
      return;
    }
  int ntotal;
  MPI_Get_count(&status,MPI_DOUBLE,&ntotal);
  vector<double> temp(ntotal);
  MPI_Recv(&temp[0],ntotal,MPI_DOUBLE,dest,tag,comm,&status);
  int n=ntotal/gradlength;
  vec.resize(n);
  for(int i=0;i<n;++i)
    {
      vec[i].density.x=temp[gradlength*i];
      vec[i].density.y=temp[gradlength*i+1];
      vec[i].pressure.x=temp[gradlength*i+2];
      vec[i].pressure.y=temp[gradlength*i+3];
      vec[i].xvelocity.x=temp[gradlength*i+4];
      vec[i].xvelocity.y=temp[gradlength*i+5];
      vec[i].yvelocity.x=temp[gradlength*i+6];
      vec[i].yvelocity.y=temp[gradlength*i+7];
      for(int j=0;j<(gradlength-8)/2;++j)
	{
	  Vector2D vtemp(temp[gradlength*i+j+8],temp[gradlength*i+j+9]);
	  vec[i].tracers.push_back(vtemp);
	}
    }
}

void SendRecvVelocity(vector<Vector2D> &vel,vector<vector<int> > sentcells,
		      vector<int> sentprocs,int totalpoints)
{
  vector<Vector2D> toadd;
  // Take care of self send velocity
  int nlist=(int)sentprocs.size();
  for(int i=0;i<nlist;++i)
    {
      if(sentprocs[i]==-1)
	{
	  if(sentcells[i].empty())
	    continue;
	  vector<Vector2D> temp=VectorValues(vel,sentcells[i]);
	  toadd.insert(toadd.end(),temp.begin(),temp.end());
	}
    }

  // Talk with other procs
  int rank,worldsize;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&worldsize);
  const vector<int> procorder=GetProcOrder(rank,worldsize);

  int n=worldsize-1;
  vector<int> sentprocs2(sentprocs);
  vector<int> talkedto;
  vector<vector<Vector2D> > vmany(sentprocs.size());
  for(int k=0;k<3;++k)
    {
      for(int i=0;i<n;++i)
	{
	  nlist=(int)sentprocs.size();
	  int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
	  // Do we talk with this processor?
	  if(index<nlist)
	    {
	      if(rank<procorder[i])
		{
		  vector<Vector2D> temp=VectorValues(vel,sentcells[index]);
		  MPI_VectorSend_Vector2D(temp,procorder[i],0,MPI_COMM_WORLD);
		  MPI_VectorRecv_Vector2D(temp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD);
		  if(!temp.empty())
		    vmany[talkedto.size()]=temp;
		}
	      else
		{
		  vector<Vector2D> temp;
		  MPI_VectorRecv_Vector2D(temp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD);
		  if(!temp.empty())
		    vmany[talkedto.size()]=temp;
		  temp=VectorValues(vel,sentcells[index]);
		  MPI_VectorSend_Vector2D(temp,procorder[i],0,MPI_COMM_WORLD);
		}
	      sentcells.erase(sentcells.begin()+index);
	      sentprocs.erase(sentprocs.begin()+index);
	      talkedto.push_back(procorder[i]);
	    }
	}
    }
  if(!vmany.empty())
    {
      int nn=(int)sentprocs2.size();
      for(int i=0;i<nn;++i)
	{
	  //int index=Find(talkedto.begin(),talkedto.end(),sentprocs2[i])-talkedto.begin();
	  if(!vmany[i].empty())
	    toadd.insert(toadd.end(),vmany[i].begin(),vmany[i].end());
	  //talkedto.erase(talkedto.begin()+index);
	  //vmany.erase(vmany.begin()+index);
	}
    }
  if(!toadd.empty())
    {
      vel.resize((size_t)totalpoints-toadd.size());
      vel.insert(vel.end(),toadd.begin(),toadd.end());
    }
}

void SendRecvGrad(vector<ReducedPrimitiveGradient2D> &grads,
		  vector<vector<int> > sentcells,vector<int> sentprocs,
		  int totalpoints)
{
  if(grads.empty())
    return;
  int gradlength=2*(int)grads[0].tracers.size()+8;
  vector<ReducedPrimitiveGradient2D> toadd;
  // Take care of self send grads
  int nlist=(int)sentprocs.size();
  for(int i=0;i<nlist;++i)
    {
      if(sentprocs[i]==-1)
	{
	  if(sentcells[i].empty())
	    continue;
	  vector<ReducedPrimitiveGradient2D> temp=VectorValues(grads,sentcells[i]);
	  toadd.insert(toadd.end(),temp.begin(),temp.end());
	}
    }

  // Talk with other procs
  int rank,worldsize;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&worldsize);
  const vector<int> procorder=GetProcOrder(rank,worldsize);

  int n=worldsize-1;
  vector<vector<ReducedPrimitiveGradient2D> > pmany(sentprocs.size());
  vector<int> sentprocs2(sentprocs);
  vector<int> talkedto;
  for(int k=0;k<3;++k)
    {
      for(int i=0;i<n;++i)
	{
	  nlist=(int)sentprocs.size();
	  int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
	  // Do we talk with this processor?
	  if(index<nlist)
	    {
	      if(rank<procorder[i])
		{
		  MPI_SendVectorGrad(VectorValues(grads,sentcells[index]),procorder[i],
				     0,MPI_COMM_WORLD);
		  vector<ReducedPrimitiveGradient2D> temp;
		  MPI_RecvVectorGrad(temp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,
				     gradlength);
		  if(!temp.empty())
		    pmany[talkedto.size()]=temp;
		}
	      else
		{
		  vector<ReducedPrimitiveGradient2D> temp;
		  MPI_RecvVectorGrad(temp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,
				     gradlength);
		  if(!temp.empty())
		    pmany[talkedto.size()]=temp;
		  MPI_SendVectorGrad(VectorValues(grads,sentcells[index]),procorder[i],
				     0,MPI_COMM_WORLD);
		}
	      sentprocs.erase(sentprocs.begin()+index);
	      sentcells.erase(sentcells.begin()+index);
	      talkedto.push_back(procorder[i]);
	    }
	}
    }
  if(!pmany.empty())
    {
      int nn=(int)sentprocs2.size();
      for(int i=0;i<nn;++i)
	{
	  //	int index=Find(talkedto.begin(),talkedto.end(),sentprocs2[i])-talkedto.begin();
	  if(!pmany[i].empty())
	    toadd.insert(toadd.end(),pmany[i].begin(),pmany[i].end());
	  //talkedto.erase(talkedto.begin()+index);
	  //pmany.erase(pmany.begin()+index);
	}
    }
  if(!toadd.empty())
    {
      grads.resize((size_t)totalpoints-toadd.size());
      grads.insert(grads.end(),toadd.begin(),toadd.end());
    }
}


int MPI_SendVectorConserved(vector<Conserved> const& vec,int dest,int tag,
			    MPI_Comm comm)
{
  if(vec.empty())
    {
      double temp=0;
      return MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
    }
  int n=(int)vec.size();
  vector<double> tosend(n*4);
  for(int i=0;i<n;++i)
    {
      tosend[4*i]=vec[i].Mass;
      tosend[4*i+1]=vec[i].Energy;
      tosend[4*i+2]=vec[i].Momentum.x;
      tosend[4*i+3]=vec[i].Momentum.y;
    }
  return MPI_Send(&tosend[0],4*n,MPI_DOUBLE,dest,tag,comm);
}

void MPI_RecvVectorConserved(vector<Conserved> &vec,int dest,int tag,
			     MPI_Comm comm)
{
  vec.clear();
  MPI_Status status;
  MPI_Probe(dest,tag,comm,&status);
  if(tag==1)
    {
      double temp;
      MPI_Recv(&temp,1,MPI_DOUBLE,dest,tag,comm,&status);
      return;
    }
  int nrecv;
  MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
  vector<double> temp(nrecv);
  MPI_Recv(&temp[0],nrecv,MPI_DOUBLE,dest,tag,comm,&status);
  int ntotal=nrecv/4;
  vec.reserve(ntotal);
  for(int i=0;i<ntotal;++i)
    {
      Conserved ctemp(temp[4*i],Vector2D(temp[4*i+2],temp[4*i+3]),temp[4*i+1]);
      vec.push_back(ctemp);
    }
}

#endif

void PeriodicUpdateCells(vector<Primitive> &cells,vector<vector<double> > &tracers,
			 vector<size_t> &customevolutions,vector<vector<int> > const& sentcells,
			 int npoints)
{
  if(sentcells.empty())
    return;
  bool traceractive=tracers.empty() ? false : true;
  if(traceractive)
    traceractive=tracers[0].empty() ? false : true;
  int n=(int)sentcells.size();
  int totalsent=0;
  for(int i=0;i<n;++i)
    if(!sentcells[i].empty())
      totalsent+=(int)sentcells[i].size();
  cells.resize(npoints-totalsent);
  customevolutions.resize(npoints-totalsent);
  if(traceractive)
    tracers.resize(npoints-totalsent);
  for(int i=0;i<n;++i)
    {
      if(sentcells[i].empty())
	continue;
      vector<Primitive> ptemp=VectorValues(cells,sentcells[i]);
      vector<size_t> ctemp=VectorValues(customevolutions,sentcells[i]);
      if(!ptemp.empty())
	cells.insert(cells.end(),ptemp.begin(),ptemp.end());
      if(!ctemp.empty())
	customevolutions.insert(customevolutions.end(),ctemp.begin(),ctemp.end());
      if(traceractive)
	{
	  vector<vector<double> > temp=VectorValues(tracers,sentcells[i]);
	  if(!temp.empty())
	    tracers.insert(tracers.end(),temp.begin(),temp.end());
	}
    }
}

void PeriodicVelocityExchange(vector<Vector2D> &vel,
			      vector<vector<int> > const& sentcells,int npoints)
{
  int n=(int)sentcells.size();
  int totalsent=0;
  for(int i=0;i<n;++i)
    if(!sentcells[i].empty())
      totalsent+=(int)sentcells[i].size();
  vel.resize(npoints-totalsent);

  for(int i=0;i<n;++i)
    {
      if(sentcells[i].empty())
	continue;
      vector<Vector2D> temp=VectorValues(vel,sentcells[i]);
      if(!temp.empty())
	vel.insert(vel.end(),temp.begin(),temp.end());
    }

}

void PeriodicGradExchange(vector<ReducedPrimitiveGradient2D> &grad,
			  vector<vector<int> > const& sentcells,int npoints)
{
  int n=(int)sentcells.size();
  int totalsent=0;
  for(int i=0;i<n;++i)
    if(!sentcells[i].empty())
      totalsent+=(int)sentcells[i].size();
  grad.resize(npoints-totalsent);

  for(int i=0;i<n;++i)
    {
      if(sentcells[i].empty())
	continue;
      vector<ReducedPrimitiveGradient2D> ptemp=VectorValues(grad,sentcells[i]);
      if(!ptemp.empty())
	grad.insert(grad.end(),ptemp.begin(),ptemp.end());
    }
}
