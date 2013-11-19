
#include <vector>

#include "comp_llpmatrix.h"
#include "lapackUtil.h"
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

struct LessdPoint {
public:

    LessdPoint() {
    };

    bool operator()(const dPoint& p1, const dPoint & p2) {
        for (int i = 0; i < p1.dimension(); i++) {
            if (p1[i] < p2[i]) {
                return true;
            } else if ((p1[i] > p2[i])) {
                return false;
            }
        }
        return false;
    }
};


extern "C" void dspev_(char &jobz, char &uplo,
        long int &n, double *ap,
        double *w, double *z,
        long int &ldz, double *work, long int &info);

void est_tangent_space_lapack(const dPoint& pt, const vector<dPoint>& neighbor_pts, double h, unsigned int tdim, vector<dVector>& tspace) {
    unsigned int dim = pt.dimension();
    assert(tdim <= dim);

    double* mat = new double[(dim + 1) * dim / 2];
    assert(mat != NULL);
    double* wr = new double[dim];
    assert(mat != NULL);
    double* vr = new double[dim * dim];
    assert(vr != NULL);

    double hh = h * h;
    memset(mat, 0, sizeof (double) * (dim + 1) * dim / 2);
    for (unsigned int i = 0; i < neighbor_pts.size(); i++) {
        dVector vec = neighbor_pts[i] - pt;
        double wt = exp(-CGAL::to_double(vec * vec) / hh);
        for (unsigned int j = 0; j < dim; j++) {
            for (unsigned int k = j; k < dim; k++) {
                mat[(k + 1) * k / 2 + j] += wt * vec[j] * vec[k];
            }
        }
    }

    // find out the eigenvalues and eigen vectors.
    char jobz = 'V', uplo = 'U';
    long int n = dim, ldz = dim, info = 0;
    double *work = new double[3 * dim];
    assert(work != NULL);

    dspev_(jobz, uplo, n, mat, wr, vr, ldz, work, info);

    //cout<<"dspev_ info: "<<info<<endl;

    for (unsigned int i = 1; i <= tdim; i++) {
        tspace.push_back(dVector(dim, vr + (dim - i) * dim, vr + (dim - i + 1) * dim));
    }


    //cout<<"wr: ";
    //for(unsigned int j = 0; j < dim; j ++){
    //	cout<<wr[j]<<" ";
    //}
    //cout<<endl;

    delete []mat;
    delete []wr;
    delete []vr;
}

extern "C" void dgetrf_(long int &m, long int &n,
        double *a, long int &lda,
        long int* ipiv, long int &info);

double simplex_sign_volume(vector<dPoint> points) {
    assert(points.size() > 0);
    unsigned int dim = points[0].dimension();

    assert(points.size() == dim + 1);

    double *A = new double[dim * dim];
    assert(A != NULL);
    for (unsigned int i = 0; i < dim; i++) {
        dVector vec = points[i + 1] - points[0];
        for (unsigned int j = 0; j < dim; j++) {
            A[i * dim + j] = CGAL::to_double(vec[j]);
        }
    }

    long int m = dim, n = dim, lda = dim, info = 0;
    long int *ipiv = new long int[dim];
    assert(ipiv != NULL);

    dgetrf_(m, n, A, lda, ipiv, info);

    double volume = 1;
    double factor = 1;
    for (unsigned int i = 0; i < dim; i++) {
        volume *= A[i * dim + i];
        factor *= (i + 1);
    }
    volume /= factor;

    delete []ipiv;
    delete []A;

    //--------------------------------------------------
    // //debug
    // for(unsigned int i = 0; i <= dim; i ++){
    // 	cout<<points[i]<<endl;
    // }
    // cout<<"volume: "<<volume<<endl;
    // getchar();
    //--------------------------------------------------

    return volume;
}


void generate_sym_pcd_qb(
        PCloud& pcloud,
        double h,
        double rho,
        vector<unsigned int>& II,
        vector<unsigned int>& JJ,
        vector<double>& SS,
        vector<double>& BB) {
    printf("generate_pcdlaplace_symmetric_matrix_sparse\n");

    std::vector< std::vector<double> > projected_v;
    std::vector<double> projected_vv;
    projected_vv.resize(2);
    double vc_area;

    typedef CGAL::Search_traits_d<KCd> Traits;
    typedef CGAL::Fair<Traits> Fair;
    typedef CGAL::Kd_tree<Traits, Fair> Tree;
    typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
    typedef KCd::FT FT;

    const unsigned int tdim = 2;


    printf("h: %f, rho: %f\n", h, rho);
    double hh = h * h;
    double nmfactor = M_PI * hh * hh / 4;

    printf("nmfactor: %f\n", nmfactor);

    unsigned int np = pcloud.p_count();
    vector<dPoint> points;
    map<const dPoint, int, LessdPoint> pt2index;
    for (unsigned int i = 0; i < np; i++) {
        points.push_back(pcloud.point(i).coord());
        pt2index.insert(make_pair(pcloud.point(i).coord(), i));
    }
    Fair fair(10);
    Tree tree(points.begin(), points.end(), fair);
    points.clear();



    double eps = 1e-3 * h * rho;
    vector<dPoint> neighbor_pts;
    vector<int> neighbor_indices;
    vector<dVector> tspace;
    Triangulation_2 tr;

    // calculate for each point p
    for (unsigned int i = 0; i < np; i++) {

        printf("i: %d\r", i);

        dPoint pt = pcloud.point(i).coord();

        //yangliu: h seems to be \eps/\rho
        //search the points within distance h * rho;
        neighbor_pts.clear();
        Fuzzy_sphere fs(pt, h * rho, FT(eps));
        tree.search(back_insert_iterator<vector<dPoint> >(neighbor_pts), fs);

        //look up the indices of the neighboring points
        neighbor_indices.clear();

        bool find = false;
        for (
                vector<dPoint>::iterator iter = neighbor_pts.begin();
                iter != neighbor_pts.end(); iter++) {
            //cout<<*iter<<endl;
            map<const dPoint, int, LessdPoint>::iterator
            iter_pt2index = pt2index.find(*iter);
            if (iter_pt2index != pt2index.end()) {
                neighbor_indices.push_back(iter_pt2index->second);
                if (iter_pt2index->second == (int) i) {
                    find = true;
                }
            } else {
                //cout<<"Warning: Fail to find the index of the neighboring points!!"<<endl;
                printf("generate_sym_pcd_qb: Failed to find the index of the neighboring points\n");
            }
        }
        assert(neighbor_indices.size() >= 3 && find);
        if ((neighbor_indices.size() < 3) || !find) printf("neighbor_indices.size() < 3 \n");
        //estimate normal
        tspace.clear();
        est_tangent_space_lapack(pt, neighbor_pts, h, tdim, tspace);

        //cout<<"pt: "<<pt<<endl;
        //cout<<"tspace: "<<tspace[0]<<", "<<tspace[1]<<endl;

        projected_v.clear();
        projected_vv[0] = 0;
        projected_vv[1] = 0;
        projected_v.push_back(projected_vv);

        tr.clear();
        for (unsigned int j = 0; j < neighbor_pts.size(); j++) {
            dVector vec = neighbor_pts[j] - pt;

            double x = CGAL::to_double(vec * tspace[0]);
            double y = CGAL::to_double(vec * tspace[1]);

            projected_vv[0] = x;
            projected_vv[1] = y;

            Vertex_handle_2d vh = tr.insert(Point_2d(x, y));
            assert(vh != NULL);
            if (vh->get_origin_id() >= 0) {
                //cout<<"warning: containing duplicated points"<<endl;
                printf("Containing duplicated points\n");
            }
            vh->set_origin_id(neighbor_indices[j]);

            projected_v.push_back(projected_vv);
        }
        // calculate correspoinding voronoi diagram cell area
        vc_area = voronoi_cell_area(projected_v);

        int v2id;
        dPoint v2;

        /**
         * Diagonal elements is not calculated. Only gaussian kernel is calculated. 
         */
        for (
                vector<dPoint>::iterator iter = neighbor_pts.begin();
                iter != neighbor_pts.end(); iter++) {
            
            map<const dPoint, int, LessdPoint>::iterator
            iter_pt2index = pt2index.find(*iter);
            if (iter_pt2index != pt2index.end()) {
                
                v2id = iter_pt2index->second;
                v2 = iter_pt2index->first;
                if (iter_pt2index->second != (int) i) {
                    double sqdist =
                            CGAL::to_double((v2 - pt) * (v2 - pt));
                    double weight = exp(-sqdist / hh) / nmfactor;

                    II.push_back(i + 1);
                    JJ.push_back(v2id + 1);
                    SS.push_back(weight);
                }
            } else {
                printf("generate_sym_pcd_qb: v2=pt in calculation.\n");
            }
        }
        /*
        for (Triangulation_2::Vertex_iterator viter = tr.vertices_begin(); viter != tr.vertices_end(); viter++) {
            unsigned int vid = viter->get_origin_id();
            if (vid == i) {

                assert((viter->point() - Point_2d(0, 0)) * (viter->point() - Point_2d(0, 0)) < FT(1e-5));
                continue;
            }

            double sqdist =
                    CGAL::to_double((viter->point() - Point_2d(0, 0)) * (viter->point() - Point_2d(0, 0)));
            double weight = exp(-sqdist / hh) / nmfactor;



            II.push_back(i + 1);
            JJ.push_back(vid + 1);
            SS.push_back(weight);
        }
         */

        BB.push_back(vc_area);
    }
    printf("\n");
}
