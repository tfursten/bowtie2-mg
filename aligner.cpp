#include <map>
#include <set>
#include <string>
#include <tclap/CmdLine.h>
#include <fstream>

std::vector<std::string> ingest_read_ids(const std::string &db_path) {
    std::cout << "reading from " << db_path << " ..." << std::endl;

    std::vector<std::string> references;

    std::ifstream input(db_path);
    std::string line, name, content;
    while (getline(input, line).good()) {

        if (line.empty() || line[0] == '>') {
            //we've hit the next header
            if (!name.empty()) {
                //if this isn't the first header, add the sequence to the references
                references.push_back(name);
                name.clear();
            }
            if (!line.empty()) {
                name = line.substr(1);
            }
            content.clear();
        } else if (!name.empty()) {
            if (line.find(' ') != std::string::npos) { // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if (!name.empty()) { // clean up last entry that wasn't inserted
        references.push_back(name);
    }
    return references;
}

int main(int argc, const char *argv[]) {
    try {
        //setup the command line parameters
        TCLAP::CmdLine cmd("Metagenomics wrapper around bowtie2.", ' ', WRAPPER_VERSION);

        TCLAP::ValueArg<int> num_workers_arg(
                "p", "processes",
                "Number of threads to allow bowtie2 to run at a time (default: 1)", false, 1, "integer");
        cmd.add(num_workers_arg);

        TCLAP::ValueArg<std::string> result_output_arg(
                "o", "output",
                "Absolute path to desired results output (including filename)", true, "", "absolute path");
        cmd.add(result_output_arg);

        TCLAP::ValueArg<std::string> index_prefix_arg(
                "i", "index",
                "Absolute path to the index prefix", true, "", "absolute path");
        cmd.add(index_prefix_arg);

        TCLAP::ValueArg<std::string> query_path_arg(
                "q", "query", "Absolute path to query reads FASTA file.", true, "", "absolute path to FASTA");
        cmd.add(query_path_arg);

        TCLAP::ValueArg<int> edit_distance_arg(
                "e", "edits", "Number of edits (subs and indels) to tolerate in a match (default: 3).", true, 3, "number");
        cmd.add(edit_distance_arg);

        cmd.parse(argc, argv);

        std::string indexpath = index_prefix_arg.getValue();
        std::string resultspath = result_output_arg.getValue();
        std::string querypath = query_path_arg.getValue();
        int num_threads = num_workers_arg.getValue();
        int max_edit_distance = edit_distance_arg.getValue();

        std::cout << "Parsed command line arguments:\nReading database from:\t\t" << indexpath <<
        "\nWriting results to:\t\t\t"
        << resultspath << "\nAllowing " << num_threads <<
        " worker threads" << std::endl;

        std::cout << "Getting all read IDs from the query..." << std::endl;

        //get all read ids from query file
        std::map<std::string, std::set<std::pair<int, std::string>>> reads_and_tax_ids;
        std::vector<std::string> ids = ingest_read_ids(querypath);

        //insert each read id into data structure where read -> list of tax ids
        for (auto i = ids.begin(); i != ids.end(); i++) {
            reads_and_tax_ids[*i].clear(); //easy way to initialize that key in the store
        }
        ids.clear();

        std::cout << "Running alignment and parsing results..." << std::endl;

        //set up buffer and execution for bowtie2 output
        FILE *pfp;
        int status;

        char buf[BUFFER_SIZE];
        std::string output;

        std::string command("bowtie2 --quiet --no-head --no-unal --omit-sec-seq --all --end-to-end --fast -f");
        command += " --threads " + std::to_string(num_threads);
        command += " -x " + indexpath;
        command += " -U " + querypath;

        std::cout << "Running this command:\n" << command << std::endl;

        //begin bowtie2 alignment
        pfp = popen(command.c_str(), "r");

        if (pfp) {
            //declare the SAM fields for later iteration
            std::string read_id;
            std::string flag;
            std::string reference_id;
            std::string position;
            std::string mapq;
            std::string cigar;
            std::string mrnm;
            std::string mpos;
            std::string tlen;
            std::string seq;
            std::string qual;
            std::string tag;

            while (fgets(buf, BUFFER_SIZE, pfp)) {
                //capture bowtie2 alignment
                output += buf;

                //see if we can find the newline
                std::string::size_type found_pos = output.find('\n');
                while (found_pos != std::string::npos) {
                    //get the first line of the input buffer
                    std::string sam_string = output.substr(0, found_pos);
                    //clip output string to remove first line (including newline)
                    output = output.substr(found_pos + 1);

                    //get ready to get all the fields from the sam string
                    std::stringstream samstream(sam_string);

                    //parse SAM output
                    std::getline(samstream, read_id, '\t');
                    std::getline(samstream, flag, '\t');
                    std::getline(samstream, reference_id, '\t');
                    std::getline(samstream, position, '\t');
                    std::getline(samstream, mapq, '\t');
                    std::getline(samstream, cigar, '\t');
                    std::getline(samstream, mrnm, '\t');
                    std::getline(samstream, mpos, '\t');
                    std::getline(samstream, tlen, '\t');
                    std::getline(samstream, seq, '\t');
                    std::getline(samstream, qual, '\t');

                    // the tags are still tab-separated, but there are a variable amount of them
                    int edit_dist = 100; //an untenable default value just in case bowtie2 doesn't emit the tag
                    while (std::getline(samstream, tag, '\t')) {
                        if (tag[0] == 'N' && tag[1] == 'M') { // NM:i: is the edit distance tag that bowtie 2 outputs
                            std::string::size_type last_colon_pos = tag.find(':', 3);
                            edit_dist = std::stoi(tag.substr(last_colon_pos + 1));
                            break;
                        }
                    }

                    // if match quality is high enough, add that tax id to the data structure
                    if (edit_dist <= max_edit_distance) {
                        std::string taxid = reference_id.substr(reference_id.find('-') + 1);
                        reads_and_tax_ids[read_id].insert(std::pair<int, std::string>(edit_dist, taxid));
                    }

                    //find the next newline if we have one
                    found_pos = output.find('\n');
                }
            }

            //close bowtie2 pid
            status = pclose(pfp);

            if (status == 0) {
                std::cout << output << std::endl;
            } else {
                std::cout << "ERROR: Non-zero exit code (" << status << ")" << std::endl;
            }
        } else {
            //TODO throw some error at the user
        }


        //for each read id
        //write to results file with list of matching tax ids

        std::ofstream resultstream(resultspath);
        for (auto i = reads_and_tax_ids.begin(); i != reads_and_tax_ids.end(); i++) {
            if (i->second.size() > 0) {
                resultstream << i->first << ':';

                int current_edit_dist = 0;
                for (auto j = i->second.begin(); j != i->second.end(); j++) {
                    while (j->first != current_edit_dist) {
                        resultstream << ':';
                        current_edit_dist++;
                    }

                    resultstream << j->second << ',';
                }
                resultstream << '\n';
            }
        }
        resultstream.close();

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}
