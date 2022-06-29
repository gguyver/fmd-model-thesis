# Readme

### Description
This repo compiles into an epidemiological model, intended to model the spread and maintenance of Foot-and-Mouth disease (FMD) in the Republic of Turkey.

It requires data on villages/farms, and recorded cattle shipments to simulate this.

---
### Inputs

The data to be imported, and options to specify aspects of the disease and control policies implemented are specified in a YAML configuration file. By default, this configuration file is named config.yaml and will be looked for in the working directory of the executable. Other names and/or locations can be specified as a CLI options, e.g. `./model-binary ../inputs/different-config.yaml`.
Three files act as inputs to the model.
1. The config.yaml file
2. A csv file describing the nodes (epi-units/villages) to be modelled, `node_file` in config.yaml.
3. A csv file describing the shipments between nodes, `shipment_file` in config.yaml

The node file should be in the format:

| ID | km north (double) | km east (double) | number of cattle (integer) |

Shipment records should be in the format:

| Day (integer) |  source node ID | target node ID | number shipped (integer) |

---
### Outputs
Outputs are written to the directory specified by `output_id`.

By default, a csv file in the tidy format is output called `$(output_id)-node-infection-report.csv`. This file describes the incidence and prevalence for each combination of:
- run
- day
- serotype
- true/detected

Optionally other descriptive files can also be output, as described in the following table:

| Filename | Configuration option | Description |
|----------|----------------------|-------------|
| $(output_id)-infection-events.csv | `infection_events` | A list of every transmission spread between nodes, including type of transmission event. |
| $(output_id)-global-compartment-sums.csv | `compartment_sums` | A sum of the global population in each compartment, for every run, day and serotype. |
| $(serotype)-$(run)-detailed-report.csv | `detailed_reports` | Reports the time since infected for each node, for each serotype and run. Output to $(output_directory)/detailed-reports | 

---
### Configuration
Configuration of the model takes place in the configfile. The model will look for a file named `config.yaml` in the same directory as the executable.
 Optionally, specifying a configuration file of another name and/or in another directory can be done via the CLI, such as:
`./<model> <path-to-different-config>` etc, so long as the normal options are present.

Make sure to preserve the structure of the file. E.g. `node_file:` should be under `inputs:`.

| Option | Type | Description |
|--------|------|-------------|
| **`inputs:`** | - | ------------------------ |
| `node_file:` | String | Filename describing nodes to run on, as described in Inputs section. |
| `shipment_file:` | String | Filename describing shipments to run, as described in Inputs section. |
| **`outputs:`** | - | ------------------------- |
| `output_directory:` | String | Directory where all outputs will be written. |
| `output_id` | String | ID which all outputs will be prefixed with, to identify the relevant results. |
| `detailed_report:` | Boolean | Whether to output detailed reports. |
| `compartment_sums:` | Boolean | Whether to output global compartment sums. |
| `infection_events:` | Boolean | Whether to output infection events. |
| `verbosity:` | Integer | Console verbosity level 0-2, 0 is silent, 2 outputs everything. |
| **`general:`** | - | ------------------------- |
| `num_runs:` | Integer | The number of replicates of the model run. Should be > 0. |
| `start_day:` | Integer | The integer day to start each model run on. |
| `end_day:` | Integer | The integer day to end each model run on. Inclusive. |
| **`model_options:`** | - | ------------------------- |
| `run_shipments:` | Boolean | Whether user-provided animal shipments will happen. Note it will not read `shipments_file` if this is set to false. |
| `implement_burn_in:` | Boolean | Whether a burn-in run should happen. Burn-in is a discarded model run from which all normal runs will then proceed. |
| `burn_in_duration:` | Integer | If `burn_in_implemented`, how many days should the burn in period be? |
| `maternal_immunity:` | Boolean | Should maternal immunity be modelled? |
| `implement_force_infections:` | Boolean | Should there be a background force of infection which can possibly infect 1 uninfected node per day? |
| `forcing_rate:` | Double | The daily probability of a forced infection. |
| `max_node_infection_duration:` | Integer | How long can a Node realistically be infected for? Possibly will be deprecated, as this is only used when Nodes are remaining infected for completely unrealistic amounts of time. |
| **`seeding:`** | - | ------------------------- |
| `fixed_seed:` | Boolean | If true, reads from `seed_nodes`. If false, randomly seed `number_random_seeds` nodes |
| `number_random_seeds:` | Integer | The number of nodes to randomly seed if `fixed_seed` is false. Should be in interval [1, number of nodes]| 
| `seed_nodes:` | List of Strings | A list of the Node IDs which should be seeded with infection at the beginning of each run (or burn-in). |
| `serotype:` | String | The serotype that each seed node will be seeded with. This needs to be present in disease->serotypes as well. |
| **`population_parameters:`** | - | ------------------------- |
| `birth_rate:` | Double | The daily rate at which births happen. |
| `mortality_rate:` | Double | The daily rate at which mortality unrelated to infection happens. |
| **`disease:`** | - | ------------------------- |
| `serotypes:` | List of Strings | The names of the serotypes to be modelled. |
| `$(serotype)_parameters:` | - | There should be one of these sections for each serotype specified in disease->serotypes. Each of the following parameters is specific to $(serotype). |
| `beta:` | Double | The within-node rate of infection. |
| `sigma:` | Double | The within-node rate of progression to infectious. Reciprocal of serotype's average latent period. |
| `gamma:` | Double | The within-node rate of recovery. Reciprocal of average infectious period. |
| `immunity_duration:` | Double | The average duration (in days) of natural immunity, i.e. from exposure to the actual disease. |
| `maternal_immunity_duration:` | Double | The average duration (in days) of maternal immunity. |
| `mortality:` | Double | The daily rate of infection related mortality. |
| `transmission:` | Double | The per-capita force of infection between-nodes. |
| `susceptibility:` | Double | The per-capita susceptiblity to transmission. Should be deprecated. |
| `shipments:` | - | Parameters for transmission via animal shipments. |
| `fomite_transmission_probability:` | Double | The percentage probability that disease is transmitted via fomites when a shipment originates from an infected Node. E.g. 1.0%, 5.0% etc... |
| **`kernel:`** | - | ------------------------- |
| `a:` | Integer | Kernel parameter. Do not change from 1. |
| `b:` | Double | Kernel parameter (shape or size?). |
| `c:` | Double | Kernel shape or size parameter. Needs to be 2.0 or greater. |
| **`control_options:`** | - | ------------------------- |
| `detection:` | - | -------------------------|
| `detection_probability:` | Double | The per-infection probability of detecting an infection. Should be between 0 (cannot detect) and 100 (detect everything). |
| `detection_delay:` | Integer | The delay between a node becoming infected and the infection being detected, assuming it is detectable. |
| `vaccine:` | - | ------------------------- |
| `$(serotype):` | - | There should be a section describing the following parameters for each serotype in disease->serotypes. |
| `efficacy:` | Double | The average efficacy of the vaccine in producing protective immunity against this serotype. 0-100. |
| `efficacy_stdev:` | Double | The standard deviation to the efficacy. |
| `duration:` | Double | The average number of days that the protective immunity generated against $(serotype) lasts. |
| `daily_vaccination_capacity:` | Integer | The number of nodes which can be vaccinated in a day. |
| `reactive_vaccination:` | - | ------------------------- |
| `implement:` | Boolean | Should reactive vaccination be modelled? |
| `radius:` | Double | What radius around the detected infected node should be vaccinated? |
| `coverage:` | Double | What percentage of those in the radius are actually vaccinated? |
| `mass_vaccination:` | - | ------------------------- |
| `implement:` | Boolean | Should mass vaccination be modelled? |
| `interval:` | Integer | If mass vaccination is modelled, how many days are there between rounds of mass vaccination? |
| `coverage:` | Double | What percentage of nodes identified for mass vaccination are actually vaccinated? |
| `movement_ban_policy:` | - | ------------------------- |
| `implement:` | Boolean | Should a reactive movement ban around detected nodes be modelled? |
| `radius:` | Double | What radius around a detected node are banned from shipping cattle in or out? |
| `duration:` | Integer | How long does a movement ban last? |
| `compliance:` | Double | What percentage of nodes in the radius comply with the ban? |