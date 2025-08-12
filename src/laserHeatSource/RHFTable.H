// RHFTable.H  (header-only)
#pragma once
#include "List.H"
#include "scalar.H"
#include "error.H"
#include "OSspecific.H"
#include <fstream>
#include <sstream>
#include <cstdlib> // add this

namespace Foam {

class RHFTable
{
    List<scalar> t_;   // seconds
    List<scalar> v_;   // RHF (dimensionless)

public:
    void readCSV(const fileName& f)
    {
        t_.clear(); v_.clear();

        Info << "RHFTable: trying to open file: " << f << endl;
        std::ifstream in(f.c_str());
        if (!in.good())
        {
            FatalErrorInFunction << "Cannot open RHF CSV: " << f << exit(FatalError);
        }

        std::string line; bool headerSkipped = false;
        while (std::getline(in, line))
        {
            Info << "RHFTable: raw line: [" << line << "]" << endl;

            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);

            Info << "RHFTable: trimmed line: [" << line << "]" << endl;

            if (line.empty()) continue;

            if (!headerSkipped) { Info << "RHFTable: skipping header line" << endl; headerSkipped = true; continue; } // skip header

            std::stringstream ss(line);
            std::string tok;
            scalar t, r;

            std::getline(ss, tok, ','); t = ::atof(tok.c_str());
            std::getline(ss, tok, ','); /* x unused */
            std::getline(ss, tok, ','); /* y unused */
            std::getline(ss, tok, ','); r = ::atof(tok.c_str());

            Info << "RHFTable: parsed t=" << t << " r=" << r << endl;

            t_.append(t);
            v_.append(r);
        }

        Info << "RHFTable: total rows read: " << t_.size() << endl;

        if (t_.empty())
        {
            FatalErrorInFunction << "RHF CSV has no rows" << exit(FatalError);
        }
    }

    // clamp at ends, linear interpolate inside
    scalar value(const scalar tq) const
    {
        if (tq <= t_.first()) return v_.first();
        if (tq >= t_.last())  return v_.last();

        label lo = 0, hi = t_.size()-1;
        while (hi - lo > 1)
        {
            const label mid = (lo + hi)/2;
            if (tq < t_[mid]) hi = mid; else lo = mid;
        }
        const scalar f = (tq - t_[lo])/(t_[hi] - t_[lo]);
        return v_[lo] + f*(v_[hi] - v_[lo]);
    }
};

} // namespace Foam

