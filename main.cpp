#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <random>

struct grasp_struct{
    int min = 65500;
    int max = 0;
    float avg = 0.0;
};

// Remove Redundant Sets
// w - weight of each subset i
// a - element j served by subset i
// b - subset i serves element j
// x - x(i) = 1 if subset i is used, 0 otherwise
int remove_redundacy(std::vector<int> &w, std::vector<std::vector<int>> &a, std::vector<std::vector<int>> &b, std::vector<bool> &x){

    while (1){

        std::vector<int> white_list; // List os subsets that are the only serving one element
        std::vector<int> black_list; // List os subsets that are serving elements already served

        // Initializes black list with all used subsets
        for(int i=0; i<x.size(); i++)
            if(x[i])
                black_list.push_back(i);

        // Checks the white list of sets - subsets that cannot be removed
        for(int j=0; j<a.size(); j++)
            if((a[j].size() == 1) && (std::find(white_list.begin(), white_list.end(), a[j][0]) == white_list.end()))
                white_list.push_back(a[j][0]);

        // Updates black list depending on the white list
        std::sort(std::begin(white_list), std::end(white_list));
        black_list.erase(std::remove_if(std::begin(black_list), std::end(black_list),
            [&](auto x){return std::binary_search(std::begin(white_list), std::end(white_list),x);}), std::end(black_list));

        // Break if no redundant sets are left
        if(!black_list.size())
            return 0;

        // Removes Redundant Sets with the largest weight per number of elements served
        float max_weight = 0.0;
        int removed_set = -1;
        for(int s=0; s<black_list.size(); s++)
            if(w[black_list[s]]/float(b[black_list[s]].size()) > max_weight){
                max_weight = w[black_list[s]]/float(b[black_list[s]].size());
                removed_set = black_list[s];
            }

        // Removes selected set
        x[removed_set] = false;

        // Update a
        for (int j=0; j<a.size(); j++)
            a[j].erase(std::remove(a[j].begin(), a[j].end(), removed_set), a[j].end());
    }
}

// Gets the total weight of the used subsets
// x - x(i) = 1 if subset i is used, 0 otherwise
// w - weight of each subset i
int get_weight(std::vector<bool> &x, std::vector<int> &w){

    int weight = 0;
    for (int i=0; i<x.size(); i++)
        weight += x[i]*w[i];

    return weight;
}

// Does Local Search of the Neighbourhood (n = 1) of the solution
// x - x(i) = 1 if subset i is used, 0 otherwise
// w - weight of each subset i
// a - element j served by subset i
// b - subset i serves element j
// first_improvement:
//   = true - FIRST improvement Local Search
//   = false - BEST improvement Local Search
std::vector<bool> local_search(std::vector<bool> x, std::vector<int> &w,  std::vector<std::vector<int>> a, std::vector<std::vector<int>> &b){

    std::vector<bool> selected_x_neighbour = x;
    std::vector<std::vector<int>> selected_a_neighbour = a;

    // Sequence of Sets Iterator
    std::vector<int> it(x.size());
    std::iota(std::begin(it), std::end(it), 0);

    // Shuffle Iterator
    std::shuffle(std::begin(it), std::end(it), std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));

    int weight = get_weight(x, w);

    bool shake = false;

    while(1){

        if(!shake){ // Neighbourhood search (N=1)
            for(int i=0; i<x.size(); i++)if(!x[it[i]]){

                    // Neighbour of the solution
                    std::vector<bool> x_neighbour = x;
                    std::vector<std::vector<int>> a_neighbour = a;

                    // Setting not used subset to true
                    x_neighbour[it[i]] = true;

                    // Saves the subsets that serve each element
                    for(int k=0; k<b[it[i]].size(); k++)
                        a_neighbour[b[it[i]][k]].push_back(it[i]);

                    // Remove Redundant Sets
                    remove_redundacy(w, a_neighbour, b, x_neighbour);

                    // Compares Neighbour Weight
                    if(weight > get_weight(x_neighbour, w)){
                        weight = get_weight(x_neighbour, w);
                        selected_x_neighbour = x_neighbour;
                        selected_a_neighbour = a_neighbour;
                    }
            }
        }
        else { // Tries to shake the current solution by adding 2.5% of non used subsets

            int prev = 0, cnt = 0;

            for(int l=0; l<x.size()*0.025; l++){//0.025

                // Neighbour of the solution
                std::vector<bool> x_neighbour = x;
                std::vector<std::vector<int>> a_neighbour = a;

                int j;                  //0.04
                for(j=prev; (cnt<x.size()*0.04 + prev) && (j<x.size()); j++)if(!x[it[j]]){

                    // Setting not used subset to true
                    x_neighbour[it[j]] = true;

                    // Saves the subsets that serve each element
                    for(int k=0; k<b[it[j]].size(); k++)
                        a_neighbour[b[it[j]][k]].push_back(it[j]);

                    cnt++;
                }
                prev = j;

                // Remove Redundant Sets
                remove_redundacy(w, a_neighbour, b, x_neighbour);

                // Compares Neighbour Weight
                if(weight > get_weight(x_neighbour, w)){
                    weight = get_weight(x_neighbour, w);
                    selected_x_neighbour = x_neighbour;
                    selected_a_neighbour = a_neighbour;

                    shake = false;

                    //break;
                }

                // Break if it already tried all possible solutions
                if(j >= x.size())
                    break;
            }
        }

        if(x == selected_x_neighbour){

            if(shake)
                break;
            else
                shake = true;
        }
        else{
            // Updates Solution
            x = selected_x_neighbour;
            a = selected_a_neighbour;
        }
    }

    return selected_x_neighbour;
}

int main(){

    std::string dir = "/home/rafa/HM/SCP/SCP-Instances/";

    // Get files name
    std::ifstream scp_files (dir + "scp.txt");
    if (!scp_files.is_open())
        std::cout << "ERROR reading scp file!" << std::endl;

    std::vector<std::string> scp; // Every instaces of scp problem (file name)
    std::string instance; // each instance of scp problem (file name)
    std::vector<int> opt; // Every instaces optimum value
    int opt_val; // optimum value

    // Save each instance of scp problem (file name)
    while(scp_files >> instance >> opt_val){
        scp.push_back(instance);
        opt.push_back(opt_val);
    }

    // Seed for Random number generator
    srand(time(nullptr));

    // Report file
    std::ofstream report(dir + "report.csv");
    std::ofstream exe_time(dir + "exe_time.csv");

    report << "inst.\tOPT.\tCH0\tBI0\tCH1\tBI1-\tBI1+\tBI1~\t" << std::endl;
    exe_time << "inst.\tCH0\tBI0\tCH1\tBI1\t" << std::endl;

    // Coefficient that represents the size of the best neighbourhood
    float alpha = 0.002; //0.02|0.002

    // Number of GRASP iterations
    int grasp_it = 500;

    auto begin = std::chrono::high_resolution_clock::now();

    for(int f=0; f<scp.size(); f++){

        report << scp[f] << "\t" << opt[f] << "\t";
        exe_time << scp[f] << "\t";

        // Read File
        std::ifstream myfile (dir + scp[f] + ".txt");
        if (!myfile.is_open())
            std::cout << "ERROR reading file!" << std::endl;

        int n; // Number of elements (j)
        int m; // Number of subsets (i)

        // Get n - number of elements
        // Get m - number of subsets
        myfile >> n >> m;

        std::vector<int> w(m); // weight of subsets i

        // Get w - weights of subsets
        for(int i=0; i<m; i++)
            myfile >> w[i];

        std::vector<std::vector<int>> a(n); // element j served by subset i
        std::vector<std::vector<int>> b(m); // subset i serves element j

        for(int j=0; j<n; j++){

            int subset_size;
            myfile >> subset_size;
            a[j].resize(subset_size);

            // Get a - subsets j..subset_size that serve element i
            // Get b
            for(int i=0; i<subset_size; i++){
                myfile >> a[j][i];
                b[a[j][i]-1].push_back(j);
            }
        }

        myfile.close();

        std::vector<std::vector<int>> a_backup = a; // Copys a to a_backup
        std::vector<std::vector<int>> b_backup = b; // Copys b to b_backup

        // Max, Min and AVG weight for CH3 and CH4 after GRASP
        grasp_struct grasp_weight;

        // Initialize chrono for GRASP
        auto chBI_t = std::chrono::microseconds::zero();

        for(int ch=0; ch<=grasp_it; ch++){

            auto start = std::chrono::high_resolution_clock::now();
            auto init = std::chrono::high_resolution_clock::now();

            // Reset a
            a.clear();
            a.resize(n);

            b = b_backup; // Copys b_backup to b

            // Representation of the Solution
            std::vector<bool> x(m, false); // x(i) = 1 if subset i is used, 0 otherwise
            std::vector<bool> y(n, false); // y(j) = 1 if element j is covered, 0 otherwise
            int n_elements = 0; // Number of elements already covered

            int t = 0;
            while(n_elements < n){

                t++;

                int used_subset = -1; // Subset selected at each iteration

                // CH1: Selecting the subset that contains the largest number of uncovered elements
                if(!ch){
                    float max_size = 0.0;
                    for(int i=0; i<m; i++)
                        if (float(b[i].size())/w[i] > max_size){
                            max_size = float(b[i].size())/w[i];
                            used_subset = i;
                        }
                }
                // CH2: Randomly selects one of the top subsets that contains the largest number of uncovered elements
                else{
                    int neigbourhood_size = std::max(int((0.2/float(t*t*t*t) + alpha) * m), 1);
                    std::vector<int> top(neigbourhood_size);

                    std::vector<float> subset_coeff(m);
                    for(int i=0; i<m; i++)
                        subset_coeff[i] = float(b[i].size())/w[i];

                    for(int j=0; j<neigbourhood_size; j++){
                        auto it = max_element(std::begin(subset_coeff), std::end(subset_coeff));
                        top[j] = std::distance(std::begin(subset_coeff), it);
                        subset_coeff[top[j]] = -1.0;
                    }
                    used_subset = top[rand() % neigbourhood_size];
                }

                // Update the number of elements already covered
                n_elements += b[used_subset].size();

                // Saves the subsets that serve each element
                for(int k=0; k<b_backup[used_subset].size(); k++)
                    a[b_backup[used_subset][k]].push_back(used_subset);

                // Setting used subset to true
                x[used_subset] = true;

                // Update the elements already covered
                for(int k=0; k<b[used_subset].size(); k++)
                    y[b[used_subset][k]] = true;

                // Erase the already covered elements from the subset list
                for(int i=0; i<m; i++)
                    for(int j=0; j<b[i].size(); j++)
                        for(int k=0; k<b[used_subset].size(); k++)
                            if((b[i][j] == b[used_subset][k]) && !((i == used_subset) && (k == j))){
                                b[i].erase(b[i].begin() + j);
                                j--;
                                break;
                            }

                // Erase the elements from the used subset
                while(b[used_subset].size())
                    b[used_subset].erase(b[used_subset].begin());
            }

            auto ch_t = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
            start = std::chrono::high_resolution_clock::now();

            // Remove Redundant Sets
            remove_redundacy(w, a, b_backup, x);
            auto chRE_t = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);

            if(ch <= 1){
                // Gets the total weight of the used subsets (After Redundancy Elimination)
                report << get_weight(x, w) << "\t";
                exe_time << (ch_t.count() + chRE_t.count())/1000.0 << "\t";
                start = std::chrono::high_resolution_clock::now();
            }

            if(!ch){
                // Local Search BEST improvement (After Redundancy Elimination)
                // Gets the total weight of the used subsets After Local Search (After Redundancy Elimination)
                std::vector<bool> x_neigbour = local_search(x, w, a, b_backup);
                report << get_weight(x_neigbour, w) << "\t";
                auto chBI_t = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
                exe_time << (ch_t.count() + chBI_t.count())/1000.0 << "\t";
                start = std::chrono::high_resolution_clock::now();
            }
            else{
                // Local Search BEST improvement (After Redundancy Elimination)
                // Gets the total weight of the used subsets After Local Search (After Redundancy Elimination)
                std::vector<bool> x_neighbour = local_search(x, w, a, b_backup);
                int weight = get_weight(x_neighbour, w);
                grasp_weight.avg += weight;
                grasp_weight.max = std::max(grasp_weight.max, weight);
                grasp_weight.min = std::min(grasp_weight.min, weight);

                chBI_t += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);

                if(ch == grasp_it){
                    grasp_weight.avg /= grasp_it;
                    report << grasp_weight.min << "\t";
                    report << grasp_weight.max << "\t";
                    report << grasp_weight.avg << "\t";
                    exe_time << (ch_t.count() + chRE_t.count() + chBI_t.count())/1000.0 << "\t";

                    // Initialize grasp_weight values
                    grasp_weight.min = 65500;
                    grasp_weight.max = 0;
                    grasp_weight.avg = 0.0;
                }
            }
            std::cout << "Instance " << scp[f] << "(CH" << ch << ") Finished. Execution Time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - init).count()/1000000.0 << " s" << std::endl;
        }
        report << std::endl;
        exe_time << std::endl;
    }

    std::cout << "Total Completion Time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - begin).count()/1000000.0 << " s" << std::endl;

    report.close();
    exe_time.close();
    return 0;
}
