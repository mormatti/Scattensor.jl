### Some rules

- Each struct should have always the attributes with a definite type (when an attribute has a type Any is not so nice... because it is difficult to initialize in general)
- Each struct should have attributes which possibly changes for every objects (for instance, no sense to make 1000 objects with the same attribute...)
- Each struct should have attributes which are preferrably not structured types (otherwise the risk is to fall to the crazy syndrome of creating large tree objects) but this is not really strict
- Each file should contain all the dispatch of a single function, and the file name should be the one of the function.
