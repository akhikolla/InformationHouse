# gretel 0.0.1

```gretel``` is a package for generalized path analysis. It implements new measures that
are described in a forthcoming research publication, as well as some traditional path 
value measures described in *Yang, Knoke* (2001).

The exported functions in this package are

* ```gpv```, which calculates Generalized Path Value
* ```ppv```, which calculates Probabilistic Path Value
* ```binary_distance```, ```peay_path_value```, ```flament_path_length```,
  ```peay_average_path_value```, and ```flament_average_path_length```, which
  calculate path value measures described in *Yang, Knoke* (2001).
* ```generate_proximities```, which generates a matrix of values representing the
  measures of optimal paths from each source node (row index) to each target node
  (column index).
* ```opt_gpv```, which identifies the path of optimal Generalized Path Value from
  a particular source node to a particular target node
* ```opt_ppv```, which identifies the path of optimal Probabilistic Path Value from
  a particular source node to a particular target node
* ```all_opt_gpv```, which identifies the 'gpv'-optimal paths from every source node
  to every target node
* ```all_opt_ppv```, which identifies the 'ppv'-optimal paths from every source node
  to every target node
* ```unpack```, which unpacks the Dijkstra-format encoded shortest paths returned by
  ```all_opt_gpv``` and ```all_opt_ppv```. See their help pages for details.

The package also includes three example datasets

* ```YangKnoke01```
* ```OpsahlEtAl10```
* ```BuchDarrah19```