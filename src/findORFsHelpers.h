#ifndef PKG_findORFs_H
#define PKG_findORFs_H

std::vector<int> orfs_as_vector(
    std::string& main_string,
    std::string s, std::string e,
    bool longestORF,
    int minimumLength);

#endif
