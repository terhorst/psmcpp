#ifndef STITCHABLE_H
#define STITCHABLE_H

class Stitchable
{
    public:
        Stitchable(double* const rho_vals, const std::vector<int> stitch_to_block) :
            rho_vals(rho_vals), stitch_to_block(stitch_to_block) {}

    protected:
    double* const rho_vals;
    const std::vector<int> stitch_to_block;
    double* map_to_rho(const int i)
    {
        size_t bk = std::distance(stitch_to_block.begin(), 
                std::upper_bound(stitch_to_block.begin(), stitch_to_block.end(), i));
        return rho_vals + bk - 1;
    }
};

#endif
