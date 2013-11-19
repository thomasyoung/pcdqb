#ifndef __PCDWRAPPER_H__
#define __PCDWRAPPER_H__

#include <vector>

class pcdwrapper {
public:
    pcdwrapper(void);
    ~pcdwrapper(void);

public:
    bool readMFile(const char *filename);
    bool readMatFile(const char *filename);

    void calculate_pcd();

    std::vector<std::vector<double> > m_points;
    std::vector<unsigned int> m_II, m_JJ;
    std::vector<double> m_SS;
    std::vector<double> m_val;

    std::vector<std::vector<unsigned int> > m_sorted_piv;

    std::vector<double> multiplyMatrix();
};
#endif
