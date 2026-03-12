#include <iostream>
#include <fstream>
#include <cstdint>
#include <sstream>
#include <string>

using namespace std;

bool extract_id(const std::string& line, uint64_t& id)
{
    if (line.empty())
        return false;

    // isolate first column
    size_t tab = line.find('\t');
    const std::string field = (tab == std::string::npos) ? line : line.substr(0, tab);

    // find suffix after '.'
    size_t dot = field.find('.');
    if (dot == std::string::npos)
        return false;

    const char* start = field.c_str() + dot + 1;

    char* endptr = nullptr;
    id = strtoull(start, &endptr, 10);

    if (endptr == start)
        return false;

    return true;
}

size_t count_lines_minus_header(const string& path)
{
    ifstream f(path);
    if (!f) {
        cerr << "Cannot open " << path << endl;
        exit(1);
    }

    size_t count = 0;
    string line;

    while (getline(f, line))
        count++;

    if (count > 0) count--; // remove header
    return count;
}

size_t intersect_count(const string& fileA, const string& fileB)
{
    ifstream A(fileA);
    ifstream B(fileB);

    if (!A || !B) {
        cerr << "Cannot open files:\n" << fileA << "\n" << fileB << endl;
        exit(1);
    }

    string lineA, lineB;
    uint64_t idA = 0, idB = 0;

    bool hasA = false;
    bool hasB = false;

    while (getline(A, lineA)) {
        if (extract_id(lineA, idA)) {
            hasA = true;
            break;
        }
    }

    while (getline(B, lineB)) {
        if (extract_id(lineB, idB)) {
            hasB = true;
            break;
        }
    }

    size_t matches = 0;

    while (hasA && hasB)
    {
        if (idA == idB) {
            matches++;

            hasA = false;
            hasB = false;

            while (getline(A, lineA)) {
                if (extract_id(lineA, idA)) {
                    hasA = true;
                    break;
                }
            }

            while (getline(B, lineB)) {
                if (extract_id(lineB, idB)) {
                    hasB = true;
                    break;
                }
            }
        }
        else if (idA < idB) {
            hasA = false;
            while (getline(A, lineA)) {
                if (extract_id(lineA, idA)) {
                    hasA = true;
                    break;
                }
            }
        }
        else {
            hasB = false;
            while (getline(B, lineB)) {
                if (extract_id(lineB, idB)) {
                    hasB = true;
                    break;
                }
            }
        }
    }

    return matches;
}

int main(int argc, char* argv[])
{
    if (argc != 6) {
        cerr << "Usage:\n";
        cerr << "contacts_stats <raw_dir> <filtered_dir> <rna_dir> <dna_dir> <output.tsv>\n";
        return 1;
    }

    string raw_dir = argv[1];
    string filt_dir = argv[2];
    string rna_dir = argv[3];
    string dna_dir = argv[4];
    string out_file = argv[5];

    string raw_rna_unique =
        raw_dir + "/" + rna_dir + "/raw_contacts_Unique_RNA.tab.rc";

    string raw_dna_unique =
        raw_dir + "/" + dna_dir + "/raw_contacts_Unique_RNA.tab.rc";

    string raw_dna_other =
        raw_dir + "/" + dna_dir + "/raw_contacts_Other.tab.rc";

    string filt_rna_unique =
        filt_dir + "/" + rna_dir + "/filtered_raw_contacts_Unique_RNA.tab.rc";

    string filt_dna_unique =
        filt_dir + "/" + dna_dir + "/filtered_raw_contacts_Unique_RNA.tab.rc";

    /* ---- Word counts ---- */

    size_t OY = count_lines_minus_header(raw_rna_unique);
    size_t OX = count_lines_minus_header(raw_dna_unique);

    /* ---- Raw contacts joins ---- */

    size_t raw_UU = intersect_count(raw_rna_unique, raw_dna_unique);
    size_t raw_UM = intersect_count(raw_rna_unique, raw_dna_other);

    /* ---- Filtered contacts joins ---- */

    size_t filt_UU = intersect_count(filt_rna_unique, filt_dna_unique);
    size_t filt_UM = intersect_count(filt_rna_unique, raw_dna_other);

    /* ---- Output TSV ---- */

    ofstream out(out_file);

    out << "category\trna_dir\tdna_dir\tU_rna/U_rna-U_dna\tU_dna/U_rna-M_dna\n";

    out << "Uniq mapped\t"
        << rna_dir << "\t"
        << dna_dir << "\t"
        << OY << "\t"
        << OX << "\n";

    out << "Raw contacts\t"
        << rna_dir << "\t"
        << dna_dir << "\t"
        << raw_UU << "\t"
        << raw_UM << "\n";

    out << "Filtered contacts\t"
        << rna_dir << "\t"
        << dna_dir << "\t"
        << filt_UU << "\t"
        << filt_UM << "\n";

    return 0;
}
