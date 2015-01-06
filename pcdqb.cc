/*
 * File:   pcdqb.cc
 * Author: yxl072100
 *
 * Created on February 6, 2009, 3:58 PM
 */

#include <stdlib.h>
#include <stdio.h>

#include "lapackUtil.h"
#include "comp_llpmatrix.h"

#include "pcdwrapper.h"

#include <string>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include "datastructure.h"

#define bufsize 1024 * 1024 * 8
/*
 *
 */

int readMatrix(char *filename, std::vector<std::vector<double>> &matrix) {
  FILE *fp;

  if ((fp = fopen((const char *)filename, "r")) == NULL) {
    printf("readAreaList: Open Input File ERROR!\n");
    return -1;
  }
  matrix.reserve(10000);

  char *buf = static_cast<char *>(malloc(bufsize)), *pch;
  while (fgets(buf, bufsize, fp)) {
    std::vector<double> nrow;
    nrow.clear();
    pch = strtok(buf, " \n");
    while (pch != NULL) {
      double val;
      if (sscanf(pch, "%lf", &val) > 0) {
        nrow.push_back(val);
      }

      pch = strtok(NULL, " \n");
    }
    if (nrow.size() > 0) {
      matrix.push_back(nrow);
    }
  }

  free(buf);
  fclose(fp);

  return 0;
}

void cal_sym_pcd_qb(std::vector<std::vector<double>> &m_points,
                    std::vector<unsigned int> &m_II,
                    std::vector<unsigned int> &m_JJ, std::vector<double> &m_SS,
                    std::vector<double> &m_BB, double searchfactor) {

  if (m_points.size() > 10) {
    double x, y, z;

    PCloud *pc;
    int i, size;
    double *ar;
    double avgs;

    m_II.resize(0);
    m_JJ.resize(0);
    m_SS.resize(0);
    m_BB.resize(0);

    size = m_points.size();
    ar = new double[3 * size];

    i = 0;

    std::vector<double> pt;

    BOOST_FOREACH (pt, m_points) {
      x = pt[0];
      y = pt[1];
      z = pt[2];

      ar[i] = x;
      ar[i + size] = y;
      ar[i + size * 2] = z;

      i++;
    }
    pc = new PCloud(ar, size, 3);
    avgs = pc->average_size(10);

    cerr << "cal_sym_pcd_qb to reserve for IJSB .\n";
    m_II.reserve(m_points.size() * 200);
    cerr << "reservation for I done " << m_points.size() * 180 << endl;
    m_JJ.reserve(m_points.size() * 200);
    cerr << "reservation for J done " << m_points.size() * 180 << endl;
    m_SS.reserve(m_points.size() * 200);
    cerr << "reservation for S done " << m_points.size() * 180 << endl;
    m_BB.reserve(m_points.size());
    cerr << "reservation for B done " << m_points.size() << endl;
    cerr << "cal_sym_pcd_qb reservation done.\n";
    /*
    void generate_sym_pcd_qb(
            PCloud& pcloud,
            double h,
            double rho,
            vector<unsigned int>& II,
            vector<unsigned int>& JJ,
            vector<double>& SS,
            vector<double>& BB) {
     */
    generate_sym_pcd_qb(*pc, avgs * 2 * searchfactor, 3, m_II, m_JJ, m_SS,
                        m_BB);

    /* put average area for all boundry vertices */
    double val, sum = 0;
    i = 0;

    BOOST_FOREACH (val, m_BB) {
      if (val > 0) {
        sum += val;
        ++i;
      } else
        printf("cal_sym_pcd_qb: boundry detected.\n");
    }
    sum /= i;

    for (i = 0; i < m_BB.size(); ++i) {
      if (m_BB[i] <= 0)
        m_BB[i] = sum;
    }

    /* add area factor for all non-diag elements
     * currently assume that
     * PCDLap code generate symmetric connectivity
     */
    std::vector<double> diag;
    diag.resize(m_points.size(), 0);
    for (i = 0; i < m_II.size(); ++i) {
      if (m_II[i] == m_JJ[i]) {
        cerr << "cal_sym_pcd_qb: errornous diagonal element \n";
        continue;
      }

      if (m_JJ[i] > m_II[i]) {
        m_SS[i] *= m_BB[m_JJ[i] - 1];
        m_SS[i] *= m_BB[m_II[i] - 1];
      } else {
        m_SS[i] *= m_BB[m_II[i] - 1];
        m_SS[i] *= m_BB[m_JJ[i] - 1];
      }
      diag[m_II[i] - 1] += m_SS[i];
    }
    /* add the diagonal elements */
    for (i = 0; i < diag.size(); ++i) {
      m_II.push_back(i + 1);
      m_JJ.push_back(i + 1);
      m_SS.push_back(0 - diag[i]);
    }

    delete pc;
    delete[] ar;
  }
}

/**
 * Usage:
 * pcdqb.exe model Q B
 *
 */
int main(int argc, char **argv) {

  if (argc < 4) {
    printf("not enough parameters: %d\n", argc);
    printf("Usage: pcdqb.exe model(input) QFile(output) BFile(output)\n");
    return 0;
  }
  /*
      pcdwrapper pw;
      pw.readMFile(argv[1]);

      printf("%d vertices read.\n", pw.m_points.size());
      if (pw.m_points.size() <= 0) {
          printf("Reading Error\n");
      }
  */
  double searchfactor = 1.0f;
  if (argc > 4) {
    int retval = sscanf(argv[4], "%lf", &searchfactor);
    //        printf("retval: %d %s\n", retval, argv[4]);
    //        printf("SearchFactor: %f\n", searchfactor);
  }
  if (searchfactor < 0) {
    searchfactor = 1.0f;
  }
  printf("SearchFactor: %f\n", searchfactor);

  std::vector<std::vector<double>> vertexCoords;
  vertexCoords.reserve(10000);
  readMatrix(argv[1], vertexCoords);
  printf("%lu vertices read.\n", vertexCoords.size());
  if (vertexCoords.size() <= 0) {
    printf("Reading Error\n");
    return (EXIT_FAILURE);
  }

  std::vector<unsigned int> I, J;
  std::vector<double> S, B;

  I.clear();
  J.clear();
  S.clear();
  B.clear();

  //    cal_sym_pcd_qb(pw.m_points, I, J, S, B);
  cal_sym_pcd_qb(vertexCoords, I, J, S, B, searchfactor);

  printf("saving q matrix to: %s\n", argv[2]);
  save_IJS(argv[2], I, J, S) ? printf("save sym IJS true\n")
                             : printf("save sym IJS false\n");
  printf("saving b matrix to: %s\n", argv[3]);
  writeDoubleArray(B, argv[3]);

  int np;
  //    np = pw.m_points.size();
  np = vertexCoords.size();

  double *p = NULL, *d = NULL;

  return (EXIT_SUCCESS);
}
