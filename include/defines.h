#ifndef DEFINES_H
#define DEFINES_H

#define DECLARE(type, var) type var;
#define DECLPTR(type, var) type* var;

#define ZERO(type, var) var = 0;
#define INVALID(type, var) var = -1;
#define MALLOC(type, var, size) var = new type[size];

#define CLEAR(type, var) var.clear();

#define READ(type, var, tree)                                               \
    tree->SetBranchStatus(#var, 1);                                         \
    tree->SetBranchAddress(#var, &var);
#define RARR(type, var, tree)                                               \
    tree->SetBranchStatus(#var, 1);                                         \
    tree->SetBranchAddress(#var, var);

#define CREATE(type, var, tree) tree->Branch(#var, &var);

#define VARCOPY(type, var, tree)                                            \
    var = tree->var;

#define VECCOPY(type, var, tree)                                            \
    if (tree->var != nullptr)                                               \
        std::copy(tree->var->begin(), tree->var->end(),                     \
                  std::back_inserter(var));

#define ARRCOPY(type, var, tree, size)                                      \
    if (tree->var != nullptr)                                               \
        std::copy(tree->var, tree->var + size, std::back_inserter(var));

#endif /* DEFINES_H */
