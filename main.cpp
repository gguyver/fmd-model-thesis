#include <iostream>
#include <chrono>

#include "Config.h"
#include "Grid.h"


int main(int argc, char *argv[]) {
    // arguments to take are only the config file, which defines everything else
    auto start_t = std::chrono::high_resolution_clock::now();
    // CONFIG - variables
    std::string configFileAddress{"config.yaml"};
    if (argc == 1) {
        configFileAddress = "config.yaml";
    } else if (argc == 2) {
        configFileAddress = argv[1];
    } else {
        std::cerr << "Too many arguments passed. Model takes only one argument, the config filepath.\n";
        exit(1);
    }
    const Config modelConfig(configFileAddress);
    if (modelConfig.m_verbosity == 2) {
        std::cout << "Number of arguments received: " << argc << "\nArguments received:\n";
        for (int i = 0; i < argc; ++i) {
            std::cout << argv[i] << '\n';
        }
        std::flush(std::cout);
    }

    // Base Grid which each run will work off
    Grid baseGrid(&modelConfig);

    // BURN IN
    if (modelConfig.m_burnIn.implement) {
        baseGrid.run("run_burn_in", true);
    }

    // RUNS
    for (size_t run = 1; run <= modelConfig.m_numModelRuns; ++run) {
        std::string runName{"run_" + std::to_string((run))};
        Grid runGrid = baseGrid;
        runGrid.run(runName, false);
        if (modelConfig.m_extraOutputs.nodeSaveState && run == modelConfig.m_numModelRuns) {
            runGrid.saveNodeState();
        }
    }

    // REPORTING - WRITING
    modelConfig.writeReports();

    auto end_t = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_t{end_t - start_t};
    if (modelConfig.m_verbosity >= 1) {
        std::cout << "\nProgram duration (s): " << elapsed_t.count() << "\nAverage Run Duration (s): "
                  << elapsed_t.count() / modelConfig.m_numModelRuns << '\n';
    }

    return 0;
}
//