#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <tclap/CmdLine.h>
#include <fstream>
#include <sys/stat.h>

std::unordered_map<std::string, std::string> parse_reference_file(const std::string &db_path) {
    std::cout << "reading from " << db_path << " ..." << std::endl;

    std::unordered_map<std::string, std::string> references;


    return references;
}

inline std::ifstream::pos_type get_file_size(const std::string &path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) {
        return 0;
    }

    return st.st_size;
}

int main(int argc, const char *argv[]) {
    try {
        //setup the command line parameters
        TCLAP::CmdLine cmd("Metagenomics wrapper around bowtie2 index builder.", ' ', WRAPPER_VERSION);

        TCLAP::ValueArg<std::string> work_dir_path_arg(
                "w", "workdir",
                "Absolute path to an EMPTY working directory", true, "", "absolute path");
        cmd.add(work_dir_path_arg);

        TCLAP::ValueArg<std::string> index_prefix_path_arg(
                "i", "index",
                "Absolute path to the index prefix", true, "", "absolute path");
        cmd.add(index_prefix_path_arg);

        TCLAP::ValueArg<std::string> tax_id_list_path_arg(
                "t", "taxids",
                "Absolute path to the newline separated list of tax IDs to index.", true, "", "absolute path");
        cmd.add(tax_id_list_path_arg);

        TCLAP::ValueArg<unsigned long> max_gb_per_chunk_arg(
                "g", "maxchunkgigs",
                "Number of gigabytes to limit a reference file chunk to (default 3).", false, 3, "number of GBs");
        cmd.add(max_gb_per_chunk_arg);

        TCLAP::ValueArg<std::string> source_path_arg(
                "s", "sourcefile",
                "Absolute path to FASTA file to index.", true, "", "absolute path");
        cmd.add(source_path_arg);

        cmd.parse(argc, argv);

        std::string index_prefix_path = index_prefix_path_arg.getValue();
        std::string work_dir_path = work_dir_path_arg.getValue();
        std::string tax_id_list_path = tax_id_list_path_arg.getValue();
        std::string source_path = source_path_arg.getValue();
        unsigned long max_bytes_per_chunk = max_gb_per_chunk_arg.getValue() * 1024 * 1024 * 1024;

        std::cout << "Reading from " << source_path << "\nGenerating temporary files in " << work_dir_path
        << "\nAccepting tax ids from " << tax_id_list_path << "\nWriting indices to " << index_prefix_path
        << "\nChunking into files max " << max_bytes_per_chunk << " bytes in length" << std::endl;

        std::unordered_set<std::string> desired_tax_ids;
        std::ifstream taxidfile(tax_id_list_path);
        std::string current_tax_id;
        while (std::getline(taxidfile, current_tax_id).good()) {
            desired_tax_ids.insert(current_tax_id);
        }

        int chunk_num = 0;
        std::string current_chunk_path = work_dir_path + "/" + std::to_string(chunk_num) + ".fasta";
        std::ofstream current_chunk_stream(current_chunk_path);
        std::cout << "Writing taxid sequences to " << current_chunk_path << std::endl;

        std::ifstream input(source_path);
        std::string line, header, sequence;
        while (getline(input, line).good()) {

            if (line.empty() || line[0] == '>') {
                //we've hit the next header
                if (!header.empty()) {
                    //if this isn't the first header, add the sequence to the references
                    std::string taxid = header.substr(header.find('-') + 1);

                    if (desired_tax_ids.find(taxid) != desired_tax_ids.end()) {
                        current_chunk_stream << '>' << header << '\n' << sequence << std::endl;

                        if (get_file_size(current_chunk_path) > max_bytes_per_chunk) {
                            current_chunk_stream.close();
                            chunk_num++;
                            current_chunk_path = work_dir_path + "/" + std::to_string(chunk_num) + ".fasta";
                            current_chunk_stream.open(current_chunk_path);
                            std::cout << "Writing taxid sequences to " << current_chunk_path << std::endl;
                        }
                    }

                    header.clear();
                    sequence.clear();
                }
                if (!line.empty()) {
                    header = line.substr(1);
                }
                sequence.clear();
            } else if (!header.empty()) {
                if (line.find(' ') != std::string::npos) { // Invalid sequence--no spaces allowed
                    header.clear();
                    sequence.clear();
                } else {
                    sequence += line;
                }
            }
        }
        if (!header.empty()) { // clean up last entry that wasn't inserted
            std::string taxid = header.substr(header.find('-') + 1);
            if (desired_tax_ids.find(taxid) != desired_tax_ids.end()) {
                current_chunk_stream << '>' << header << '\n' << sequence << std::endl;
            }
        }
        current_chunk_stream.close();


        for (auto i = 0; i <= chunk_num; i++) {
            std::string current_chunk = work_dir_path + "/" + std::to_string(i) + ".fasta";
            std::string current_index_prefix(index_prefix_path + std::to_string(i));
            std::string bt2index_command("bowtie2-build -f -o 3 -q " + current_chunk_path + " " + current_index_prefix);

            //run the bowtie2-build command
            FILE *fp;
            int status;
            char buf[BUFFER_SIZE];
            std::string bowtie2build_output;
            fp = popen(bt2index_command.c_str(), "r");

            while (fgets(buf, BUFFER_SIZE, fp)) {
                bowtie2build_output += buf;
            }

            status = pclose(fp);

            if (status == 0) {
                std::cout << "Successfully built index chunk " << i << " of " << chunk_num << std::endl;
            } else {
                std::cout << "ERROR: Non-zero exit code from bowtie2-build (" << status << ")!\n" << bowtie2build_output << std::endl;
            }
        }

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}
