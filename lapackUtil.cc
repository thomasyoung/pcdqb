#include "lapackUtil.h"
#include "point_cloud.h"
#include "comp_llpmatrix.h"

#include "datastructure.h"
#include <math.h>
//ARPACK++ headers


#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>
#include <string>

lapackUtil::lapackUtil(void) {
}

lapackUtil::~lapackUtil(void) {
}

// IJS starts from 1

int convertIJS2colMatrix(
        std::vector<unsigned int> II,
        std::vector<unsigned int> JJ,
        std::vector<double> SS,
        double* a, int m, int n) {
    if ((m * n <= 0) ||
            (II.size() != JJ.size()) ||
            (JJ.size() != SS.size()) ||
            (a == NULL)
            )
        return -1;

    int i, j, pos;
    j = II.size();

    for (i = 0; i < j; i++) {
        pos = (JJ[i] - 1) * m + II[i] - 1;
        a[pos] = SS[i];
    }

    return 0;
}

int convertVector2Array(std::vector<double> v, double* y) {
    double val;
    int i;

    i = 0;

    BOOST_FOREACH(val, v) {
        y[i] = val;
        ++i;
    }

    return 0;
}

std::vector<double> readDoubleArray(const char* filename) {
    std::vector<double> result;

    std::ifstream ifs;
    std::string line;

    std::string s;
    std::vector<std::string> stvec;

    double val;

    ifs.open(filename);

    if (ifs.is_open()) {
        std::getline(ifs, line);
        if (line.length() > 0) {
            boost::algorithm::split(stvec, line, boost::algorithm::is_any_of(" \n"));

            BOOST_FOREACH(s, stvec) {
                try {
                    val = boost::lexical_cast<double>(s);
                    result.push_back(val);
                } catch (boost::bad_lexical_cast&) {
                    printf("lexical_cast fail\n");
                }
            }
        }
    }

    return result;
}

bool writeDoubleArray(std::vector<double>& param, const char* filename) {
    std::ofstream ofs;
    char buf[300];
    double val;
    std::string s;

    ofs.open(filename);

    BOOST_FOREACH(val, param) {
        sprintf(buf, "%4.20f ", val);
        s.clear();
        s.append(buf);

        ofs << s;
    }

    ofs.close();

    return true;
}

double voronoi_cell_area(std::vector<std::vector<double> > &m_points_2d) {
    //cerr<<"voronoi_cell_area begin\n";

    std::vector<double> pt;
    // double val;

    VD vd;
    Point_2* pt2;
    Site_2* t;

    BOOST_FOREACH(pt, m_points_2d) {
        t = new Site_2(pt[0], pt[1]);
        vd.insert(*t);
        delete t;
    }

    if (vd.is_valid()) {
        //        printf("voronoi_cell_area: valid \n");
    } else {
        cerr << "voronoi_cell_area: not valid \n";
        return 0;
    }

    //    printf("voronoi_cell_area: number_of_vertices: %d\n", vd.number_of_vertices());

    Face_handle *fi;
    pt2 = new Point_2(m_points_2d[0][0], m_points_2d[0][1]);
    Locate_result lr = vd.locate(*pt2);
    fi = boost::get<Face_handle > (&lr);
    delete pt2;

    Ccb_halfedge_circulator ec_start = (*fi)->outer_ccb();
    Ccb_halfedge_circulator ec = ec_start;

    Polygon pl;

    CGAL::Cartesian<double>::Point_2 *p2;

    double x, y;
    do {
        if (!ec->has_source()) {
            cerr << "voronoi_cell_area : error : no source.\n";
//            delete fi;
            return 0;
        }

        x = ec->source()->point().x();
        y = ec->source()->point().y();
        //        printf("polygon vertices: %f, %f\n", x, y);

        //cerr<<"just before p2 = new CGAL::Cartesian<double>::Point_2(x, y);";
        p2 = new CGAL::Cartesian<double>::Point_2(x, y);
        pl.push_back(*p2);
        delete p2;
    } while (++ec != ec_start);

    //delete fi;

    //cerr<<"voronoi_cell_area end\n";
    return pl.area();
}

int readPurePointsFile(const char* filename, std::vector<std::vector<double> > &m_points) {
    std::ifstream ifs;
    std::string line;

    std::string s;
    std::vector<std::string> stvec;

    ifs.open(filename);

    double val;
    std::vector<double> pt;


    if (ifs.is_open()) {
        m_points.clear();
        do {
            pt.clear();

            std::getline(ifs, line);
            if (line.size() <= 0) continue;
            boost::algorithm::split(stvec, line, boost::algorithm::is_any_of(" \n"));

            BOOST_FOREACH(s, stvec) {
                try {
                    val = boost::lexical_cast<double>(s);
                    pt.push_back(val);
                } catch (boost::bad_lexical_cast&) {
                    printf("readPointsFile: lexical_cast fail at No. %d\n",
                            m_points.size());
                }
            }
            m_points.push_back(pt);
        } while (!ifs.eof());

        return 1;
    } else {
        return 0;
    }
}





bool save_IJS(
        const char* filename,
        std::vector<unsigned int>& II,
        std::vector<unsigned int>& JJ,
        std::vector<double>& SS) {

    if (filename == NULL) return false;
    if ((II.size() < SS.size()) ||
            (JJ.size() < SS.size())) {
        cerr<<"save_IJS: size of parameters do not match: "<< II.size()<< JJ.size()<< SS.size()<<endl;
        return false;
    }

    cerr<<"save_IJS: parameters valid"<<endl;
    
    std::ofstream ofs;
    char buf[300];
    double val;
    std::string s;

    ofs.open(filename);

    int i = 0;

    cerr<<"save_IJS: begin iteration"<<endl;
    BOOST_FOREACH(val, SS) {
        sprintf(buf, "%d %d %4.20f\n", II[i], JJ[i], val);
        s.clear();
        s.append(buf);

        ofs << s;
        ofs.flush();
        ++i;
    }

    cerr<<"save_IJS: end iteration"<<endl;
    
    ofs.close();

    return true;
}


/**
 * intends to encapsulate lapack routines for eigen-solving:
 GEBAL    balance matrix
 int dgebal_(char *job, integer *n, doublereal *a, integer *
        lda, integer *ilo, integer *ihi, doublereal *scale, integer *info)

 GEHRD    reduce to Hessenberg form A=QHQ^H
 int dgehrd_(integer *n, integer *ilo, integer *ihi,
        doublereal *a, integer *lda, doublereal *tau, doublereal *work,
        integer *lwork, integer *info)

 ORGHR    generate matrix Q
 int dorghr_(integer *n, integer *ilo, integer *ihi,
        doublereal *a, integer *lda, doublereal *tau, doublereal *work,
        integer *lwork, integer *info)

 HSEQR    find eigenvalues and schur factorization (QR)
 int dhseqr_(char *job, char *compz, integer *n, integer *ilo,
         integer *ihi, doublereal *h__, integer *ldh, doublereal *wr,
        doublereal *wi, doublereal *z__, integer *ldz, doublereal *work,
        integer *lwork, integer *info)

 TREVC    find eigenvectors from schur factorization
 int dtrevc_(char *side, char *howmny, logical *select,
        integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
        ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m,
        doublereal *work, integer *info)

 GEBAK    transform eigenvectors of balanced matrix to those of the original matrix
 int dgebak_(char *job, char *side, integer *n, integer *ilo,
        integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *
        ldv, integer *info)
 */

extern "C" {
    /* Subroutine */ int dgebal_(char *job, long int *n, double *a, long int *
            lda, long int *ilo, long int *ihi, double *scale, long int *info);
    /* Subroutine */ int dgehrd_(long int *n, long int *ilo, long int *ihi,
            double *a, long int *lda, double *tau, double *work,
            long int *lwork, long int *info);
    /* Subroutine */ int dorghr_(long int *n, long int *ilo, long int *ihi,
            double *a, long int *lda, double *tau, double *work,
            long int *lwork, long int *info);
    /* Subroutine */ int dhseqr_(char *job, char *compz, long int *n, long int *ilo,
            long int *ihi, double *h__, long int *ldh, double *wr,
            double *wi, double *z__, long int *ldz, double *work,
            long int *lwork, long int *info);
    /* Subroutine */ int dtrevc_(char *side, char *howmny, long int *select,
            long int *n, double *t, long int *ldt, double *vl, long int *
            ldvl, double *vr, long int *ldvr, long int *mm, long int *m,
            double *work, long int *info);
    /* Subroutine */ int dgebak_(char *job, char *side, long int *n, long int *ilo,
            long int *ihi, double *scale, long int *m, double *v, long int *
            ldv, long int *info);

    void dsygst_(long int* itype, char* uplo, long int* n, double* a,
            long int* lda, double* b, long int* ldb, long int* info);
    /* 
     DSYGST reduces a real symmetric-definite generalized eigenproblem
     to standard form.
     */
}

extern "C" {
    void dgees_(char *jobvs, char *sort, int* select, int *n,
            double *a, int *lda, int *sdim, double *wr,
            double *wi, double *vs, int *ldvs, double *work,
            int *lwork, int *bwork, int *info);
    // a simple driver that computes all or part of the Schur factorization of A, with optional ordering of the eigenvalues;
    int dgeev_(char *jobvl, char *jobvr, long int *n, double *a, long int *lda, double *wr, double *wi, double *vl,
            long int *ldvl, double *vr, long int *ldvr, double *work,
            long int *lwork, long int *info);
    // a simple driver that computes all the eigenvalues of A, and (optionally) the right or left eigenvectors (or both);
}

/**
 * Input matrix is taken as column-major matrix
 *
 */
int solve_eigen(double* v,
        double* lambda, long int n) {
    if ((NULL == v) || (NULL == lambda) || (0 >= n)) {
        return -1;
    }

    // begin calculation here
    long int info = 0;

    char jobvl = 'N', jobvr = 'V';

    long int lwork = 6 * n;
    double *work = new double[lwork];

    double *wi = new double[n];
    double *wr = new double[n];
    double *vl = new double[n];
    double *vr = new double[n * n];

    long int ldvl = 1, ldvr = n;
    /*
     *    int dgeev_(char *jobvl, char *jobvr, long int *n, double *a, long int *lda, double *wr, double *wi, double *vl,
            long int *ldvl, double *vr, long int *ldvr, double *work,
            long int *lwork, long int *info);
     *int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *
        a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl,
        integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work,
        integer *lwork, integer *info);
     */
    printf("solve_eigen: %c %c %ld %ld\n", jobvl, jobvr, n, lwork);
    dgeev_(&jobvl, &jobvr, &n, v, &n, wr, wi, vl,
            &ldvl, vr, &ldvr, work,
            &lwork, &info);
    if (0 != info) {
        printf("solve_eigen error: dgeev_ : %ld\n", info);
    }

    int piv;
    for (int i = 0; i < n; ++i) {
        lambda[i] = wr[i];
        if (wi[i] != 0) {
            printf("solve_eigen: image part found: %d, %f\n", i, wi[i]);
            // try to replace imaginary part with real part
            piv = n*i;
            memcpy(v+piv+n, v+piv, sizeof(double)*n);
        } else {
            //nothing to do. normal
        }
    }
    memcpy(v, vr, sizeof (double) * n * n);

    delete[] vr;
    delete[] vl;
    delete[] wr;
    delete[] wi;

    delete[] work;

    return 0;
}

int convert_IJS_2_col_sq_matrix(
        double* result,
        std::vector<unsigned int> II,
        std::vector<unsigned int> JJ,
        std::vector<double> SS, int n) {
    if (n <= 0) return -1;

    if ((II.size() > JJ.size()) || (II.size() > SS.size()) || (NULL == result))
        return -1;

    unsigned int row, i = 0, col;
    int piv;

    BOOST_FOREACH(row, II) {
        col = JJ[i];
        piv = (row - 1)+(col - 1) * n;
        if (piv >= 0) {
            result[piv] = SS[i];
        } else {
            printf("convert_IJS_2_col_matri error: piv<0 at line %d\n", i);
        }

        ++i;
    }
    return 0;
}

int writeDoubleRowSqMatrix(double* p, int n, const char* filename) {
    if ((NULL == p) || (NULL == filename) || (n <= 0)) return -1;

    std::ofstream ofs;
    char buf[300];
    std::string s;

    ofs.open(filename);

    int piv;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            piv = i * n + j;
            sprintf(buf, "%4.20f ", p[piv]);
            s.clear();
            s.append(buf);

            ofs << s;
        }
        sprintf(buf, "\n");
        s.clear();
        s.append(buf);

        ofs << s;
    }

    ofs.close();

    return 0;
}

int writeDoubleRow(double* p, int n, const char* filename) {
    if ((NULL == p) || (NULL == filename) || (n <= 0)) return -1;

    std::ofstream ofs;
    char buf[300];
    std::string s;

    ofs.open(filename);

    for (int i = 0; i < n; ++i) {
        sprintf(buf, "%4.20f ", p[i]);
        s.clear();
        s.append(buf);

        ofs << s;
    }

    ofs.close();

    return 0;
}

