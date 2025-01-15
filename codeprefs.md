A new struct should:

- have always the attributes with a definite type (when an attribute has a type Any is not so nice... because it is difficult to initialize in general)
- have attributes which possibly changes for every objects (for instance, no sense to make 1000 objects with the same attribute...)
- have attributes which are preferrably not structured types (otherwise the risk is to fall to the crazy syndrome of creating large tree objects) but this is not really strict
