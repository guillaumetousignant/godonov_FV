#include "entities/Face_t.h"

FVM::Entities::Face_t::Face_t() {}

FVM::Entities::Face_t::Face_t(size_t node_0, size_t node_1, size_t cell_0, size_t cell_1) : 
    nodes_{node_0, node_1},
    cells_{cell_0, cell_1} {}

FVM::Entities::Face_t::~Face_t() {}
