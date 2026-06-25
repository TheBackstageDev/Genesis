#include "atom_storage.hpp"

namespace sim
{
    atomStorage::atomStorage(size_t n) 
    {
        resize(n);
    }

    void atomStorage::resize(size_t n) 
    {
        N = n;
        x.resize(n);
        y.resize(n);
        z.resize(n);

        vx.resize(n);
        vy.resize(n);
        vz.resize(n);

        fx.resize(n);
        fy.resize(n);
        fz.resize(n);

        invMass.resize(n);
        q.resize(n);
        ljParams.resize(n * 2);

        De.resize(n);
        re.resize(n);
        a.resize(n);
    }

    void atomStorage::clear()
    {
        N = 0;

        x.clear();
        y.clear();
        z.clear();

        vx.clear();
        vy.clear();
        vz.clear();

        fx.clear();
        fy.clear();
        fz.clear();

        invMass.clear();
        q.clear();
        ljParams.clear();

        De.clear();
        re.clear();
        a.clear();
    }

    void atomStorage::clearForces()
    {
        std::fill(fx.begin(), fx.end(), 0.0f);
        std::fill(fy.begin(), fy.end(), 0.0f);
        std::fill(fz.begin(), fz.end(), 0.0f);
    }

    void atomStorage::setAtom(const atomView& aView) 
    {
        size_t i = aView.i;

        if (i >= N) return;

        x[i] = aView.px; y[i] = aView.py; z[i] = aView.pz;
        vx[i] = aView.vx; vy[i] = aView.vy; vz[i] = aView.vz;
        fx[i] = aView.fx; fy[i] = aView.fy; fz[i] = aView.fz;
        invMass[i] = aView.invMass;
        q[i] = aView.q;
        ljParams[i*2]     = aView.sigma;
        ljParams[i*2 + 1] = aView.epsilon;
        
        De[i] = aView.De;
        re[i] = aView.re;
        a[i] = aView.a;
    }
} // namespace sim
