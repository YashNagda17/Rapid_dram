#include <bits/stdc++.h>
#include <string>
#include <future>
#include <thread>
#include <fstream>
#include <ramalutor/src/config.h>
#include <src/allfiles.h>
#include <ramalutor/src/config.h>
#include <ramalutor/src/DRAM.h>

using namespace std;

void extract_val(int mode)
{
    string targetfile = '';
    if (mode==0)
    {
        targetfile = 'seq2rapid.csv';
    }
    else if (mode==1)
    {
        targetfile = 'seq2rapid_lin.csv';
    }
    else
    {
        targetfile = 'seq2rapid_aff.csv';
    }
    a = 'ramulator/raptraces/464.h264ref.gz'
    ofstream MyFile(a);
    string line;
    float mem,time;
    int currentLine = 1;
    while (getline(MyFile, line)) {
        if (currentLine == 486) {
            time = line[26:40];
        }
        else if (currentLine == 617)
        {
            memory = line[30:44];
        }
        if (currentLine >= lineNum2)
        {
            break;
        }
    }
    MyFile.close();
    vector<double> data = {time,memory}
    std::ofstream file(filename);
    for (const auto& row : data) {
        for (auto iter = row.begin(); iter != row.end(); ++iter) {
            file << *iter;
            if (std::next(iter) != row.end()) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
}

void write_res(string filename, string res)
{
    string a =  filename + '.txt';
    ofstream MyFile(a);
    MyFile << res;
    MyFile.close();
    if (filename == 'rapid')
    {
        extract_val(0);
    }
    else if (filename == 'seq2rapid_lin')
    {
        extract_val(1);
    }
    else
    {
        extract_val(2);
    }
}

void rapid_band_affine(string file1, string file2)
{
    vector<char> a;
    vector<char> b;
    long long int b = 100;
    long long int sigmao = 100;
    long long int sigmae = 1;
    int m = x.length();
    int n = y.length();
    a = read_file(file1);
    b = read_file(file2);
    vector<int> row; ///Integrate Affine Gap Penalty Code
    long long int col = 4096*(2*b);
    long long int max_loc = 16777216;
    long long int curr_loc = 0;


    for (int i = 1; i < m; i++) {
        curr_loc = max_loc + col + 3*i;
        tmp = sigmao + (i - 1) * sigmae;
        write_dram(tmp,curr_loc);
    }

    for (int i = 1; i <= m; i++) {
        curr_loc = max_loc + i + col;
        write_dram(add_dram_val(curr_loc,sigmae), curr_loc + 4);
        for (int j = max(1,j-b/2); j <= min(n,j-b/2); j++) {
            curr_loc_a = max_loc + col + (i-1);
            curr_loc_c = col*(i-1) + (j-1);
            curr_loc_b = curr_loc_a + 1;
            curr_loc_c2 = curr_loc_c - 2;
            add_write_dram(curr_loc_a - 2,curr_loc_c, curr_loc_b);
            add_write_dram(curr_loc_a, curr_loc_b, curr_loc_c2);
        }
        for (int j = max(0,j-b/2); j <= min(n,j-b/2); ++j) {
            curr_loc_a = max_loc + col + (j-1);
            curr_loc_c2 = curr_loc_c - 2;
            curr_loc_c = col*(i-1) + (j);
            curr_loc_b = a+1;
            find_min_tmp(add(curr_loc_a,curr_loc_c), add_dram_val(curr_loc_b,sigmao), curr_loc_c2);
        }
    }
    string res = find_min_path(max_loc + m + b);
    string filename = 'res2rapid';
    write_res(filename,res);
    return;
}

void rapid_band_lin(string file1, string file2)
{
    vector<char> a;
    vector<char> b;
    long long int b = 100;
    int m = x.length();
    int n = y.length();
    int penalty = 1;
    int sigma = 1;
    a = read_file(file1);
    b = read_file(file2);
    vector<int> row;
    long long int col = 4096*(b);
    long long int max_loc = 16777216;
    long long int curr_loc = 0;

    for(int i = 1; i<=m; i ++ )
    {
        for (int j = 1;j<=n;j++)
        {
            curr_loc = col*(i-1) + (j-1);
            if (a[i-1]!=b[j-1])
            {
                write_dram(penalty, curr_loc);
            }
            else
            {
                write_dram(0, curr_loc);
            }
        }
    }
    for (int i = 1; i < m; i++) {
        curr_loc = max_loc + col + 2*i;
        write_dram(sigma*i,curr_loc);
    }
    for (int i = 1; i <= m; i++) {
        curr_loc = max_loc + i + col;
        write_dram(add_dram_val(curr_loc,sigma), curr_loc + 4);
        for (int j = max(1,j-b/2); j <= min(n,j-b/2); j++) {
            curr_loc_a = max_loc + col + (i-1);
            curr_loc_c = col*(i-1) + (j-1);
            curr_loc_b = curr_loc_a + 1;
            curr_loc_c2 = curr_loc_c - 2;
            add_write_dram(curr_loc_a - 2,curr_loc_c, curr_loc_b);
            add_write_dram(curr_loc_a, curr_loc_b, curr_loc_c2);
        }
        for (int j = max(0,j-b/2); j <= min(n,j-b/2); ++j) {
            curr_loc_a = max_loc + col + (j-1);
            curr_loc_c2 = curr_loc_c - 2;
            curr_loc_c = col*(i-1) + (j);
            curr_loc_b = a+1;
            find_min_tmp(add(curr_loc_a,curr_loc_c), add_dram_val(curr_loc_b,sigma), curr_loc_c2);
        }
    }
    string res = find_min_path(max_loc + m + b);
    string filename = 'seq2rapid_lin';
    write_res(filename,res);
    return;
}

void rapid(string file1, string file2)
{
    long long int cycles_used = 0;
    vector<char> a;
    vector<char> b;
    int m = x.length();
    int n = y.length();
    int penalty = 1;
    int sigma = 1;
    a = read_file(file1);
    b = read_file(file2);
    vector<int> row;
    long long int col = 4096*(n);
    long long int max_loc = 16777216;
    long long int curr_loc = 0;

    for(int i = 1; i<=m; i ++ )
    {
        for (int j = 1;j<=n;j++)
        {
            curr_loc = col*(i-1) + (j-1);
            if (a[i-1]!=b[j-1])
            {
                write_dram(penalty, curr_loc);
            }
            else
            {
                write_dram(0, curr_loc);
            }
        }
    }
    for (int i = 1; i < m; i++) {
        curr_loc = max_loc + col + 2*i;
        write_dram(sigma*i,curr_loc);
    }
    for (int i = 1; i <= m; i++) {
        curr_loc = max_loc + i + col;
        write_dram(add_dram_val(curr_loc,sigma), curr_loc + 4);
        for (int j = 1; j <= n; j++) {
            curr_loc_a = max_loc + col + (i-1);
            curr_loc_c = col*(i-1) + (j-1);
            curr_loc_b = curr_loc_a + 1;
            curr_loc_c2 = curr_loc_c - 2;
            add_write_dram(curr_loc_a - 2,curr_loc_c, curr_loc_b);
            add_write_dram(curr_loc_a, curr_loc_b, curr_loc_c2);
        }
        for (int j = 0; j <= n; ++j) {
            curr_loc_a = max_loc + col + (j-1);
            curr_loc_c2 = curr_loc_c - 2;
            curr_loc_c = col*(i-1) + (j);
            curr_loc_b = a+1;
            add(curr_loc_a,curr_loc_c);
            find_min_tmp(add(curr_loc_a,curr_loc_c), add_dram_val(curr_loc_b,sigma), curr_loc_c2);
        }
    }
    string res = find_min_path(max_loc + m + n);
    string filename = 'seq2rapid';
    write_res(filename,res);
    return;
}


int main()
{
    run_config();
    string file1 = 'seq1file1_small';
    string file2 = 'seq1file2_small';
    //rapid(file1,file2);
    rapid_band_lin(file1,file2);
    //rapid_band_affine(file1,file2);
    return 0;
}