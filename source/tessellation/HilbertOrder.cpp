#include "HilbertOrder.hpp"

namespace
{
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
    Context (vector<int>* par,AllLocalData const& local_data,int i):
      _params(par), _data(local_data), _i(i) {}

    vector<int>* _params;
    AllLocalData _data;
    int _i;
  };

  bool point_compare_x(Vector2D const& p1,Vector2D const& p2)
  {
    return p1.x<p2.x;
  }

  bool point_compare_y(Vector2D const& p1,Vector2D const& p2)
  {
    return p1.y<p2.y;
  }

  pair<int,int> rot(int n, pair<int,int> const& origin, pair<int,int> const& r)
  {
    pair<int,int> res = origin;
    if(r.second==0){
      if(r.first==1)
	res = pair<int,int>(n-1-res.first,n-1-res.second);
      return pair<int,int>(res.second,res.first);
    }
    else
      return res;
  }

  //convert (x,y) to d
  int xy2d (int n, int x, int y) {
    int d=0;
    pair<int,int> temp(x,y);
    for (int s=n/2; s>0; s/=2) {
      const pair<int,int> r((x&s)>0,(y&s)>0);
      d += s * s * ((3 * r.first) ^ r.second);
      temp = rot(s, temp, r);
    }
    return d;
  }
}

namespace {
  vector<int> range(int n)
  {
    vector<int> res(static_cast<size_t>(n));
    for(size_t i=0;i<res.size();++i)
      res[i] = static_cast<int>(i);
    return res;
  }

  template<class T> void zero_pad_start(vector<T>& v,int n)
  {
    const vector<int> temp(static_cast<size_t>(n));
    v.insert(v.begin(),temp.begin(),temp.end());
  }

  int calc_d(pair<double, double> const& dxdy,
	     pair<double, double> const& xminymin,
	     Vector2D const& cor,
	     int pow2)
  {
    if(dxdy.first<=0&&dxdy.second<=0)
      return 0;

    const pair<int, int> index
      (dxdy.first <=0 ? 0 : static_cast<int>((cor.x-xminymin.first)/dxdy.first),
       dxdy.second <=0 ? 0 : static_cast<int>((cor.y-xminymin.second)/dxdy.second));
    return xy2d(pow2,index.first,index.second);
  }
}

vector<int> HilbertOrder(vector<Vector2D> const& cor,int num,int innernum)
{
  if(cor.empty())
    return vector<int> ();
  const int N=num-innernum;
  vector<int> insert_order;

  insert_order.reserve(static_cast<size_t>(N));
  stack<Context*> context_list;
  stack<vector<int>*> IndexTableStack;
  const int temp_n=max(static_cast<int>(sqrt(static_cast<double>(N))),2);
  const int temp_pow2=max(static_cast<int>(pow(2.0,static_cast<int>(log(1.0*temp_n)/log(2.0)))),2);
  AllLocalData d_local_data(0,N,temp_n,temp_pow2);
  vector<int> p = range(N);
  vector<int>* param = &p;
  Context *context=new Context(param,d_local_data,0);
  context_list.push(context);
  bool expand;
  const vector<Vector2D> cortemp(cor.begin()+innernum,
				 cor.begin()+num);
  while(!context_list.empty())
    {
      context=context_list.top();
      AllLocalData local_data=context->_data;
      param=context->_params;
      int i=context->_i;
      context_list.pop();
      delete context;
      if(static_cast<int>(context_list.size())>N)
	{
	  UniversalError eo("Error in creating the hilbert order, probably two points are identical");
	  throw eo;
	}
      if(i==0)
	{
	  vector<Vector2D> p_cor(static_cast<size_t>(local_data.N));
	  for(size_t j=0;j<static_cast<size_t>(local_data.N);++j)
	    p_cor[j] = cortemp[static_cast<size_t>(param->at(j))];
	  const double xmin=(*min_element(p_cor.begin(),
					  p_cor.end(),
					  point_compare_x)).x;
	  const double xmax=(*max_element(p_cor.begin(),
					  p_cor.end(),
					  point_compare_x)).x;
	  const double ymin=(*min_element(p_cor.begin(),
					  p_cor.end(),
					  point_compare_y)).y;
	  const double ymax=(*max_element(p_cor.begin(),
					  p_cor.end(),
					  point_compare_y)).y;
	  const double dx=static_cast<double>(xmax-xmin)/(local_data.pow2-1);
	  const double dy=static_cast<double>(ymax-ymin)/(local_data.pow2-1);
	  local_data.IndexTable=new vector<int>[local_data.pow2*local_data.pow2];
	  for(int j=0;j<local_data.N;j++){
	    const int d = calc_d
	      (pair<double,double>(dx,dy),
	       pair<double,double>(xmin,ymin),
	       p_cor[static_cast<size_t>(j)],local_data.pow2);
	    local_data.IndexTable[d].push_back(param->at(static_cast<size_t>(j)));
	  }
	}
      if(i<(local_data.pow2*local_data.pow2-1))
	{
	  context = new Context(param, local_data,i+1);	// i+1 is the for loop increment
	  context_list.push(context);
	}
      if(local_data.IndexTable[i].size()<2)
	{
	  expand=false;
	  for(size_t k=0;k<local_data.IndexTable[i].size();k++)
	    insert_order.push_back(local_data.IndexTable[i][k]);
	}
      else
	{
	  expand=true;
	  int temp_n2=max(static_cast<int>(sqrt(double(local_data.IndexTable[i].size()))),2);
	  context=new Context
	    (&local_data.IndexTable[i],
	     AllLocalData
	     (0,
	      static_cast<int>(local_data.IndexTable[i].size()),
	      temp_n2,
	      max(static_cast<int>(pow(2.0,static_cast<int>(log(1.0*temp_n2)/log(2.0)))),2)), 0);
	  context_list.push(context);
	}
      if(i==(local_data.pow2*local_data.pow2-1))
	{
	  for(size_t l=0;l<IndexTableStack.size();++l)
	    {
	      delete [] IndexTableStack.top();
	      IndexTableStack.pop();
	    }
	  if(!expand)
	    delete [] local_data.IndexTable;
	  else
	    IndexTableStack.push(local_data.IndexTable);
	}
    }
  if(innernum>0)
    {
      zero_pad_start(insert_order,innernum);
      for(int i=0;i<num;++i)
	{
	  if(i<innernum)
	    insert_order[static_cast<size_t>(i)]=i;
	  else
	    insert_order[static_cast<size_t>(i)]+=innernum;
	}
    }
  return insert_order;
}
