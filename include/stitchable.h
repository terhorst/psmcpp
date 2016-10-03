#ifndef STITCHABLE_H
#define STITCHABLE_H

class Stitchable
{
    public:
        Stitchable(double* const rho_vals, const std::vector<int> stitch_to_block) :
            rho_vals(rho_vals), stitch_to_block(stitch_to_block) {
                // std::cout << "RHO #1: " << rho_vals << std::endl;
                // std::cout << "RHO #2: " << rho_vals + 1 << std::endl;
                // std::cout << "RHO #3: " << rho_vals + 2 << std::endl;
                // std::cout << "RHO #1 VAL: " << *(rho_vals) << std::endl;
                // std::cout << "RHO #2 VAL: " << *(rho_vals + 1) << std::endl;
                // std::cout << "RHO #3 VAL: " << *(rho_vals + 2) << std::endl;

            }

    protected:
    double* const rho_vals;
    const std::vector<int> stitch_to_block;
    double* map_to_rho(const int i)
    {
        size_t bk = std::distance(stitch_to_block.begin(), 
                std::upper_bound(stitch_to_block.begin(), stitch_to_block.end(), i));
        if(stitch_to_block.size() == bk){
            // TODO: Stupid thing because sometimes we don't have equal length obs?? Figure out after conference
            //std::cout << "EYE VALUE" << i << std::endl;
            return rho_vals + bk - 2;
        }
        return rho_vals + bk - 1;
    }
};

#endif
