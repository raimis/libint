/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_src_lib_libint_atom_h_
#define _libint2_src_lib_libint_atom_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
# error "libint2/atom.h requires C++11 support"
#endif

#include <array>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include <tuple>

#include <libint2/chemistry/elements.h>

namespace libint2 {

  struct Atom {
      int atomic_number;
      double x, y, z;
  };
}

namespace {

bool strcaseequal(const std::string& a, const std::string& b) {
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(),
                                            [](char a, char b) {return ::tolower(a) == ::tolower(b);}
                                           );
}

/// reads the list of atoms from a file in the standard or PBC-extended XYZ format
/// \sa libint2::read_dotxyz
/// \sa libint2::read_dotxyz_pbc
inline std::tuple<std::vector<libint2::Atom>,
                  std::array<std::array<double, 3>, 3>>
__libint2_read_dotxyz(std::istream& is, const bool pbc = false) {
  using libint2::Atom;

  // first line = # of atoms
  size_t natom;
  is >> natom;

  // read off the rest of first line and discard
  std::string rest_of_line;
  std::getline(is, rest_of_line);

  // second line = comment
  std::string comment;
  std::getline(is, comment);

  // rest of lines are atoms (and unit cell parameters, if pbc = true)
  const auto nlines_expected = natom + (pbc ? 3 : 0);
  std::vector<Atom> atoms(natom, Atom{0, 0.0, 0.0, 0.0});
  std::array<std::array<double, 3>, 3> unit_cell({{0.0, 0.0, 0.0}});
  for (auto line = 0, atom_index = 0; line < nlines_expected; ++line) {
    if (is.eof())
      throw std::runtime_error(std::string("libint2::read_dotxyz: expected ") +
                               std::to_string(nlines_expected) +
                               " sets of coordinates but only " +
                               std::to_string(line) + " received");

    // read line
    std::string linestr;
    std::getline(is, linestr);
    std::istringstream iss(linestr);
    // then parse ... this handles "extended" XYZ formats
    std::string element_symbol;
    double x, y, z;
    iss >> element_symbol >> x >> y >> z;

    // .xyz files report Cartesian coordinates in angstroms; convert to bohr
    const auto angstrom_to_bohr = 1 / 0.52917721092;  // 2010 CODATA value
//    const auto angstrom_to_bohr =
//        1 / 0.529177249;  // 1986 CODATA value, used by MPQC

    auto assign_atom = [angstrom_to_bohr](Atom& atom, int Z, double x, double y,
                                          double z) {
      atom.atomic_number = Z;
      atom.x = x * angstrom_to_bohr;
      atom.y = y * angstrom_to_bohr;
      atom.z = z * angstrom_to_bohr;
    };
    auto assign_xyz = [angstrom_to_bohr](std::array<double, 3>& xyz, double x,
                                         double y, double z) {
      xyz[0] = x * angstrom_to_bohr;
      xyz[1] = y * angstrom_to_bohr;
      xyz[2] = z * angstrom_to_bohr;
    };

    auto axis = -1;
    // if pbc = true, look for unit cell params
    if (pbc) {
      if (strcaseequal("AA", element_symbol))
        axis = 0;
      if (strcaseequal("BB", element_symbol))
        axis = 1;
      if (strcaseequal("CC", element_symbol))
        axis = 2;
      if (axis != -1) {
        assign_xyz(unit_cell[axis], x, y, z);
      }
    }

    // .xyz files report element labels, hence convert to atomic numbers
    if (axis == -1) {
      int Z = -1;
      using libint2::chemistry::element_info;
      for (const auto& e : element_info) {
        if (strcaseequal(e.symbol, element_symbol)) {
          Z = e.Z;
          break;
        }
      }
      if (Z == -1) {
        std::cerr << "libint2::read_dotxyz: element symbol \"" << element_symbol
                  << "\" is not recognized" << std::endl;
        throw "Did not recognize element symbol in .xyz file";
      }

      assign_atom(atoms[atom_index++], Z, x, y, z);
    }
  }

  return std::make_tuple(atoms,unit_cell);
}

}

namespace libint2 {

  /// reads the list of atoms from a file in the standard XYZ format supported
  /// by most chemistry software
  /// \note see Wikipedia XYZ format entry at http://en.wikipedia.org/wiki/XYZ_file_format
  inline std::vector<Atom> read_dotxyz(std::istream& is) {
    std::vector<Atom> atoms;
    std::tie(atoms, std::ignore) = __libint2_read_dotxyz(is);
    return atoms;
  }

  /// reads the list of atoms from a file in the PBC-extended XYZ format
  /// \note The unit cell vectors in PBC-extended XYZ file are specified as atoms with
  ///       element symbols "AA", "BB", and "CC" (N.B. the element symbols are not
  ///       case-sensitive). If a unit cell vector is omitted, it is assumed to be zero.
  ///       Omitting all three unit cell vectors is equivalent to an infinite unit
  ///       cell (no periodicity in any direction).
  /// \param is the input stream
  /// \return a tuple composed of the list of atoms and an array of 3
  ///         unit cell vectors, \c A , \c B , and \c C .
  inline auto read_dotxyz_pbc(std::istream& is)
      -> decltype(__libint2_read_dotxyz(is, true)) {
    return __libint2_read_dotxyz(is, true);
  }

  /// converts a vector of <code>Atom</code>s to a vector of point charges
  std::vector<std::pair<
      double,
      std::array<double, 3>>> inline make_point_charges(const std::
                                                            vector<
                                                                libint2::Atom>&
                                                                atoms) {
    std::vector<std::pair<double, std::array<double, 3>>> q(atoms.size());
    for (const auto& atom : atoms) {
      q.emplace_back(static_cast<double>(atom.atomic_number),
                     std::array<double, 3>{{atom.x, atom.y, atom.z}});
    }
    return q;
  }

} // namespace libint2

#endif /* _libint2_src_lib_libint_atom_h_ */
