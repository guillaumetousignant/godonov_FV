#include "entities/Cell_t.h"

FVM::Entities::Cell_t::Cell_t(int n_sides) : 
        nodes_(n_sides),
        cells_(n_cells) {}

FVM::Entities::Cell_t::~Cell_t() {}