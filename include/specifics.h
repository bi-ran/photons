#ifndef SPECIFICS_H
#define SPECIFICS_H

template <typename T>
static int64_t within_hem_failure_region(T* t, int64_t index) {
    return ((*t->phoSCEta)[index] < -1.3
        && (*t->phoSCPhi)[index] < -0.87
        && (*t->phoSCPhi)[index] > -1.57);
}

#endif /* SPECIFICS_H */
