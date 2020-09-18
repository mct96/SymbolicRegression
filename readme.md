```bash
+-------------------------------------------------------------------------------+                                     
|                              Symbolic Regression                              |
|                              -------------------                              |
|                                                                               |
|    Author: Matheus Cândido Teixeira.                                          |
|    Date: 31/08/2020                                                           |
|    Last Update: 12/08/2020                                                    |
|    Cuiabá - MT                                                                |
|                                                                               |
+-------------------------------------------------------------------------------+
|                                     USAGE                                     |
|                                                                               |
|./sreg DATASET [-h] [-o filename.csv] [-n 1..] [--population-size 1..*]        |
|               [--max-generation 1..*] [--generation-method method]            |
|               [--max-depth 1..*] [--selection-method method] [-k 1..*]        |
|               [-pm 0.01 .. *] [-pc 0.01 .. *] [--error-metric metric]         |
|               [--fitness-threshold 0.01 .. *] [-E] --seeds [filename].csv     |
+---------------------------+--------------------------------------+------------+
|   Parameter               |   Arguments                          | Needed     |
+---------------------------+--------------------------------------+------------+
|	[input dataset].csv     | Position argument MUST be the first  |     *      |
|	-h, --help              |                                      |            |
|	-o, --output            | [filename].csv                       |            |
|	-n, --variables         | Integer                              |            |
|	--population-size       | Integer                              |            |
|	--max-generation        | Integer                              |            |
|	--generation-method     | ramped-hh, full, grow                |            |
|	--max-depth             | Integer                              |            |
|	--selection-method      | roulette, tournament                 |            |
|	-k, --k-tournament      | Integer                              |            |
|	-pm, --prob-mutation    | Real                                 |            |
|	-pc, --prob-crossover   | Real                                 |            |
|	--error-metric          | mae, mse, rmse                       |            |
|	--fitness-threshold     | Real                                 |            |
|	-E, --eletism           |                                      |            |
|	--seeds                 | [filename].csv                       |     *      |
+---------------------------+--------------------------------------+------------+
|   Repository: https://github.com/MatheusCTeixeira/SymbolicRegression          |
|   E-mail: matheuscandido2009@gmail.com, matheus.candido@dcc.ufmg.br           |
+-------------------------------------------------------------------------------+
```