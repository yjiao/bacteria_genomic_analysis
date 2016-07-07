# Visualization of genomic analysis in Matlab

Genomic analysis typically integrate many sources of data in addition to mutations. These may include drug susceptibility of the bacteria that were sequenced, as well as house keeping data such as the type of experiment and position of the strain on a 96-well plate.

These types of data typically exist as separate files and are imported into Matlab in matrix form. It is easy to manipulate data in matrix form in Matlab, but it is difficult to remember the different formats and locations of all relevant data. For example, a single strain may have many mutations, each of which are associated with different parameters that reside in different rows of different matrices. This strain may further have drug susceptibility data in a different matrix. It would be convenient for other types of analyses to bind all these data together into a single strain object.

Whether to use an object-oriented approach or to stick with matrices in Matlab depends on the type of analysis required downstream. To resolve this tension, both matrices and objects were created, and functions were written to convert one to the other.

In this case, our data can be represented by three objects:

* **Strain** contains information about a particular strain of bacteria that was sequenced. A **Strain** object keeps track of variables such as drug susceptibility, the type of experiment, and location on 96-well plates. A **Strain** can have many **Mutation** objects.
* **Mutation** contains information about a genetic variant, such as the location in the genome it is found in, the reference and alternate, etc. It also keeps track of the **Strain** objects with this particular mutation.
* **Event** objects keep track of sequencing information, such as the number of reads at the location and all breseq parameters associated with it. They represent the intersection of a **Strain** and a **Mutation**.

