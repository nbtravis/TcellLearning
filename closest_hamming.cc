#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;

const string SELF_FILE = "./data/self-9mers.txt";
const vector<string> VIRUS_FILES = {
  "./data/hepb-9mers.txt",
  "./data/hiv-9mers.txt",
  "./data/hepc-9mers.txt",
  "./data/zika-9mers.txt",
  "./data/vac-9mers.txt",
  "./data/hcmv-9mers.txt",
  "./data/lis-9mers.txt",
  "./data/mal-9mers.txt"
};

vector<string> GetPeptidesFromFile(const string& filename) {
    ifstream input_file(filename);

    vector<string> peptides;
    string peptide;
    while (getline(input_file, peptide)){
        peptides.push_back(peptide);
    }
    return peptides;
}

int HammingDist(const string& s1, const string& s2) {
  int dist = 0;
  for (int i = 0; i < s1.size(); ++i) {
    dist += static_cast<int>(s1[i] != s2[i]);
  }
  return dist;
}

int ClosestHammingToSelf(const string& peptide, const vector<string>& self_peptides) {
  int min_dist = 9;
  for (const string& self_peptide : self_peptides) {
    int dist = HammingDist(peptide, self_peptide);
    min_dist = min(min_dist, dist);
  }
  return min_dist;
}


int main() {
  vector<string> self_peptides = GetPeptidesFromFile(SELF_FILE);
  unordered_map<string, vector<string>> virus_file_to_peptides;
  for (const string& filename : VIRUS_FILES) {
    auto virus_peptides = GetPeptidesFromFile(filename);
    virus_file_to_peptides.insert({filename, virus_peptides});
  }

  for (const string& filename : VIRUS_FILES) {   
    cout << "Processing viruses in file " << filename << endl;
    
    const vector<string>& virus_peptides = virus_file_to_peptides[filename];

    int i = 0;
    int len = virus_peptides.size();
    vector<int> dist_counts = {0, 0, 0, 0, 0, 0, 0};
    for (const string& virus_peptide : virus_peptides) {
      int min_dist = ClosestHammingToSelf(virus_peptide, self_peptides);
      ++dist_counts[min_dist];
      ++i;
      if (i % 100 == 0) {
        cout << "Processing peptide " << i << " / " << len << endl;
      }
    }


    cout << "min hamming distance counts: [";
    for (int dist = 0; dist < 7; ++dist) {
      cout << dist_counts[dist];
      if (dist < 6) {
        cout << ", ";
      } else {
        cout << "]\n" << endl;
      }
    }    
  }
  
  return 0;
}  
