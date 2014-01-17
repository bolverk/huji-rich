#include "HilbertOrder.hpp"

namespace
{
  class AllMyParameters
  {
  public:

    AllMyParameters(vector<int> * vin):
      v(vin) {}

    AllMyParameters(const AllMyParameters &other):
      v(other.v) {}

    vector<int> *v;
  };

  class AllLocalData
  {
  public:

    AllLocalData(vector<int> *vin,int iN,int inn,int ipow2):
      IndexTable(vin), N(iN), n(inn), pow2(ipow2) {}

    AllLocalData(const AllLocalData &other):
      IndexTable(other.IndexTable),
      N(other.N),
      n(other.n),
      pow2(other.pow2) {}

    vector<int> *IndexTable;
    int N;
    int n;
    int pow2;
  };

  class Context
  {
  public:
    Context (const AllMyParameters &par,const AllLocalData &local_data,int i):
      _params(par), _data(local_data), _i(i) {}

    AllMyParameters _params;
    AllLocalData _data;
    int _i;
  };

  bool point_compare_x(Vector2D &p1,Vector2D &p2)
  {
    if(p1.get_x()<p2.get_x())
      return true;
    else
      return false;
  }

  bool point_compare_y(Vector2D &p1,Vector2D &p2)
  {
    if(p1.get_y()<p2.get_y())
      return true;
    else
      return false;
  }

  //rotate/flip a quadrant appropriately
  void rot(int n, int *x, int *y, int rx, int ry) {
    if (ry == 0) {
      int t;
      if (rx == 1) {
	*x = n-1 - *x;
	*y = n-1 - *y;
      }
      t  = *x;
      *x = *y;
      *y = t;
    }
  }

  //convert (x,y) to d
  int xy2d (int n, int x, int y) {
    int rx, ry, s, d=0;
    for (s=n/2; s>0; s/=2) {
      rx = (x & s) > 0;
      ry = (y & s) > 0;
      d += s * s * ((3 * rx) ^ ry);
      rot(s, &x, &y, rx, ry);
    }
    return d;
  }
}

vector<int> HilbertOrder(vector<Vector2D> const& cor,int num,int innernum)
{
  vector<Vector2D> cortemp;
  int N=num-innernum;
  cortemp.resize(N);
  copy(cor.begin()+innernum,cor.begin()+num,cortemp.begin());
  //if(innernum>0)
  //	cortemp.erase(cortemp.begin(),cortemp.begin()+innernum);
  vector<int> insert_order;
  vector<int> p;

  insert_order.reserve(N);
  p.resize(N);
  for(int i=0;i<N;++i)
    p[i]=i;
  double AvgNumInCell=1;
  int MaxOcc_upy=2;
  stack<Context*> context_list;
  vector<int> *last_IndexTable;
  stack<vector<int>*> IndexTableStack;
  int temp_n=max((int)sqrt((double)N/AvgNumInCell),2);
  int temp_pow2=max((int) pow(2.0,(int)(log(1.0*temp_n)/log(2.0))),2);
  AllLocalData d_local_data(0,N,temp_n,temp_pow2);
  AllMyParameters param(&p);
  Context *context=new Context(param,d_local_data,0);
  context_list.push(context);
  bool expand;
  while(!context_list.empty())
    {
      context=context_list.top();
      AllLocalData local_data=context->_data;
      param=context->_params;
      int i=context->_i;
      context_list.pop();
      delete context;
      if((int)context_list.size()>N)
	throw("Error in creating the hilbert order, probably two points are identical");
      if(i==0)
	{
	  vector<Vector2D> *p_cor=new vector<Vector2D>;
	  p_cor->reserve(local_data.N);
	  for(int j=0;j<local_data.N;j++)
	    p_cor->push_back(cortemp[param.v->at(j)]);
	  double xmin=(*min_element(p_cor->begin(),p_cor->end(),point_compare_x)).get_x();
	  double xmax=(*max_element(p_cor->begin(),p_cor->end(),point_compare_x)).get_x();
	  double ymin=(*min_element(p_cor->begin(),p_cor->end(),point_compare_y)).get_y();
	  double ymax=(*max_element(p_cor->begin(),p_cor->end(),point_compare_y)).get_y();
	  double dx=(double)(xmax-xmin)/(local_data.pow2-1),dy=(double)(ymax-ymin)/(local_data.pow2-1);
	  local_data.IndexTable=new vector<int>[local_data.pow2*local_data.pow2];
	  int d,xIndex,yIndex;
	  for(int j=0;j<local_data.N;j++)
	    {
	      if(dx<=0||dy<=0)
		{
		  if(dx<=0)
		    {
		      yIndex=(int)((p_cor->at(j).get_y()-ymin)/dy);
		      xIndex=0;
		    }
		  else
		    {
		      xIndex=(int)((p_cor->at(j).get_x()-xmin)/dx);
		      yIndex=0;
		    }
		}
	      else
		{
		  xIndex=(int)((p_cor->at(j).get_x()-xmin)/dx);
		  yIndex=(int)((p_cor->at(j).get_y()-ymin)/dy);
		}
	      d=xy2d(local_data.pow2,xIndex,yIndex);
	      local_data.IndexTable[d].push_back(param.v->at(j));
	    }
	  delete p_cor;
	}
      if(i<(local_data.pow2*local_data.pow2-1))
	{
	  context = new Context(param, local_data,i+1);	// i+1 is the for loop increment
	  context_list.push(context);
	}
      if((int)local_data.IndexTable[i].size()<MaxOcc_upy)
	{
	  expand=false;
	  for(int k=0;k<(int)local_data.IndexTable[i].size();k++)
	    insert_order.push_back(local_data.IndexTable[i][k]);
	}
      else
	{
	  expand=true;
	  //int temp_n2=max((int)sqrt(1.0*(local_data.IndexTable[i].size())/AvgNumInCell),2);
	  int temp_n2=max((int)sqrt(double(local_data.IndexTable[i].size())/AvgNumInCell),2);
	  int temp_pow22=max((int) pow(2.0,(int)(log(1.0*temp_n2)/log(2.0))),2);
	  AllMyParameters param2(&local_data.IndexTable[i]);
	  AllLocalData local_temp(0,int(local_data.IndexTable[i].size()),temp_n2,temp_pow22);
	  context=new Context(param2,local_temp,0);
	  context_list.push(context);
	}
      if(i==(local_data.pow2*local_data.pow2-1))
	{
	  int NN=(int)IndexTableStack.size();
	  for(int l=0;l<NN;++l)
	    {
	      last_IndexTable=IndexTableStack.top();
	      delete [] last_IndexTable;
	      IndexTableStack.pop();
	    }
	  if(!expand)
	    delete [] local_data.IndexTable;
	  else
	    {
	      last_IndexTable=local_data.IndexTable;
	      IndexTableStack.push(last_IndexTable);
	    }
	}
    }
  if(innernum>0)
    {
      vector<int> temp(innernum,0);
      insert_order.insert(insert_order.begin(),temp.begin(),temp.end());
      for(int i=0;i<num;++i)
	{
	  if(i<innernum)
	    insert_order[i]=i;
	  else
	    insert_order[i]+=innernum;
	}
    }
  return insert_order;
}
