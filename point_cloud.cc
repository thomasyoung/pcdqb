#include <math.h>
#include <queue>
#include <fstream>
//#include <GL/glut.h>

#include "point_cloud.h"

double __partcolorgrtable[31][3] = {{0, 0, 1},
                                    {0, 0.0625, 0.9375},
                                    {0, 0.1250, 0.8750},
                                    {0, 0.1875, 0.8125},
                                    {0, 0.2500, 0.7500},
                                    {0, 0.3125, 0.6875},
                                    {0, 0.3750, 0.6250},
                                    {0, 0.4375, 0.5625},
                                    {0, 0.5625, 0.4375},
                                    {0, 0.6250, 0.3750},
                                    {0, 0.6875, 0.3125},
                                    {0, 0.7500, 0.2500},
                                    {0, 0.8125, 0.1875},
                                    {0, 0.8750, 0.1250},
                                    {0, 0.9375, 0.0625},

                                    {0, 1, 0},
                                    {0.0625, 0.9375, 0},
                                    {0.1250, 0.8750, 0},
                                    {0.1875, 0.8125, 0},
                                    {0.2500, 0.7500, 0},
                                    {0.3125, 0.6875, 0},
                                    {0.3750, 0.6250, 0},
                                    {0.4375, 0.5625, 0},
                                    {0.5625, 0.4375, 0},
                                    {0.6250, 0.3750, 0},
                                    {0.6875, 0.3125, 0},
                                    {0.7500, 0.2500, 0},
                                    {0.8125, 0.1875, 0},
                                    {0.8750, 0.1250, 0},
                                    {0.9375, 0.0625, 0},
                                    {1, 0, 0}};

void get_color(double v, double minv, double maxv, double c[3]) {
  double scale = 0.8;
  /*
    //---------------------------------------------------------------------------------------
    //all warm to cold
    double nnv = maxv - minv;
    int minc = 0, maxc = 63;
    double vs = (v - minv) / nnv * (maxc - minc + 1);

    int i = (int)(vs);
    for(int j = 0; j < 3; j ++)
      c[j] = scale * ((vs - i) * (__colortable[i + 1][j] - __colortable[i][j]) +
    __colortable[i][j]);

    //  cout<<"c: "<<c[0]<<" "<<c[1]<<" "<<c[2]<<endl;
    //---------------------------------------------------------------------------------------
  */

  // maxv = 2.5;
  // minv = 1;
  double nnv = maxv - minv;

  if (fabs(nnv) < MYNZERO) {
    for (int j = 0; j < 3; j++)
      c[j] = scale * __partcolorgrtable[15][j];
    return;
  }

  //---------------------------------------------------------------------------------------
  // warm to cold
  int minc = 0, maxc = 30;
  if (v <= minv) {
    for (int j = 0; j < 3; j++)
      c[j] = scale * __partcolorgrtable[minc][j];
    return;
  }

  if (v >= maxv) {
    for (int j = 0; j < 3; j++)
      c[j] = scale * __partcolorgrtable[maxc][j];
    return;
  }
  double vs = (v - minv) / nnv * (maxc - minc);

  int i = (int)(vs);
  for (int j = 0; j < 3; j++)
    c[j] =
        scale *
        ((vs - i) * (__partcolorgrtable[i + 1][j] - __partcolorgrtable[i][j]) +
         __partcolorgrtable[i][j]);
  //  cout<<"c: "<<c[0]<<" "<<c[1]<<" "<<c[2]<<endl;
  //---------------------------------------------------------------------------------------

  /*
    //---------------------------------------------------------------------------------------
    //grey scale

    double maxc = 0, minc = 0.5;

    if(v <= minv) { c[0] = c[1] = c[2] = minc; return; }
    if(v >= maxv) { c[0] = c[1] = c[2] = maxc; return; }

    double t = (v - 1) / nnv;
    c[0] = c[1] = c[2] = minc  + (maxc - minc) * t;

    //---------------------------------------------------------------------------------------
  */
}

PCloud::PCloud(double *points, unsigned int np, unsigned int dim) {
  vector<double> coords;
  coords.resize(dim);
  for (unsigned int i = 0; i < np; i++) {
    for (unsigned int j = 0; j < dim; j++) {
      coords[j] = points[j * np + i];
    }
    add_point(VPCloud(dPoint(dim, coords.begin(), coords.end())));
  }
  // GetBBox();
}

bool PCloud::ReadPointCloud(char *filename) {
  ifstream fin;
  fin.open(filename);
  if (!fin) {
    cout << "Cannot open the file " << filename << endl;
    return false;
  }

  istream_iterator<double> input(fin);
  istream_iterator<double> end;

  if (input == end) {
    cout << "Empty data " << endl;
    return false;
  }

  unsigned int d = (unsigned int)(*input);
  ++input;
  set_dd(d);
  cout << "d: " << dd() << endl;

  vector<double> coord;
  coord.resize(d);
  while (input != end) {
    for (unsigned int i = 0; i < d - 1; i++) {
      coord[i] = *input;
      ++input;
      // cout<<coord[i]<<" ";

      if (input == end) {
        cout << "Invalid data" << endl;
        return false;
      }
    }
    coord[d - 1] = *input;
    ++input;
    // cout<<coord[d - 1]<<endl;
    add_point(VPCloud(dPoint(d, coord.begin(), coord.end())));
  }
  fin.close();
  GetBBox();
  return true;
}

////////////////////////////////////////////////////////////////////////////////////////
// GetBBox
//--------
// Get a box which is twice as bigger as the boundding box of the  object
//
void PCloud::GetBBox() {
  typedef KCd::FT FT;

  if (p_count() == 0) {
    return;
  }
  vector<FT> min;
  min.resize(dd(), FT(FLT_MAX));
  vector<FT> max;
  max.resize(dd(), FT(-FLT_MAX));

  for (unsigned int i = 0; i < p_count(); i++) {
    for (unsigned int j = 0; j < dd(); j++) {
      if ((point(i).coord())[j] < min[j]) {
        min[j] = (point(i).coord())[j];
      } else if ((point(i).coord())[j] > max[j]) {
        max[j] = (point(i).coord())[j];
      }
    } // for j
  } // for i

  _pmin = dPoint(dd(), min.begin(), min.end());
  _pmax = dPoint(dd(), max.begin(), max.end());

  cerr << "BBox: min " << _pmin << " max " << _pmax << endl;
}

double PCloud::average_size(unsigned int k) {
  typedef CGAL::Search_traits_d<KCd> Traits;
  typedef CGAL::Euclidean_distance<Traits> Distance;
  typedef CGAL::Fair<Traits> Fair;
  typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance, Fair>
      Neighbor_search;
  typedef Neighbor_search::Tree Tree;
  typedef CGAL::Fair<Traits> Fair;
  typedef KCd::FT FT;

  unsigned int np = p_count();

  vector<dPoint> points;
  for (unsigned int i = 0; i < np; i++) {
    points.push_back(point(i).coord());
  }

  Fair fair(10);
  Tree tree(points.begin(), points.end(), fair);

  double eps = 1e-3;
  double as = 0;
  for (unsigned int i = 0; i < np; i++) {
    dPoint pt = point(i).coord();
    Neighbor_search search(tree, pt, k, FT(eps));

    FT dist(0);
    unsigned int nfound = 0;
    for (Neighbor_search::iterator iter = search.begin(); iter != search.end();
         iter++) {
      dist += sqrt(iter->second);
      nfound++;
    }
    if (nfound > 0) {
      as += (dist / nfound);
    }
  }
  return as / np;
}

void PCloud::OutPCloud(char *filename) {
  ofstream fout;
  fout.open(filename);
  if (fout.fail()) {
    std::cout << "Failed to open file!" << std::endl;
    return;
  }
  for (unsigned int i = 0; i < p_count(); i++) {
    fout << point(i).coord() << endl;
  }
  fout.close();
}
