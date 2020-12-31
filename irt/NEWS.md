

# irt 0.1.3

* Added `$max_score` for `Testlet` objects. For a given `Testlet` object `t1`, 
  `t1$max_score` returns the maximum obtainable score by all items within 
  the testlet. 
* Added `$resp_item_list` for `Itempool` objects. For an `Itempool` object 
  `ip`, `ip$resp_item_list` will combine items that are not in a testlet and 
  items within a testlet and returns a list object that is composed of only
  `Item` objects. 


# irt 0.1.2

* Package can be used in MacOS as well. Resolved compiler related issues. 
* Added `est_bilog()` function to estimate item parameters for "1PL", "2PL" and
"3PL" models using BILOG-MG on computers with Windows OS. This function 
requires an installed BILOG-MG program on the machine. 


# irt 0.1.0

* This is the initial release of the package. 

