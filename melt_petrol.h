/*
    Copyright (C) 2023 Petra Maierova

    This file is a plugin for ASPECT. It is partly derived from
    Melt global material model of ASPECT.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ASPECT; see the file LICENSE.  If not, see
    <https://www.gnu.org/licenses/>.

*/

#ifndef _aspect_material_model_melt_petrol_h
#define _aspect_material_model_melt_petrol_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

#include <aspect/mesh_deformation/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that implements a simple formulation of
     * melting/freezing and different compositions of solid and melt.
     * Equilibrium melt fraction (porosity) is given by the lever rule.
     * Solidus and liquidus temperatures depend linearly on
     * (reference lithostatic) pressure and composition.
     * Model requires tracking of porosity, solid and melt (liquid) compositions.
     * Assumes cartesian coordinates and constant vertical gravity.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltPetrol : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>, public MaterialModel::MeltFractionModel<dim>
    {
    public:
      /**
       * Defines that the MeltPetrol model is incompressible.
       */
      bool is_compressible() const override;
  
      /**
       * Update the base model and viscosity function at the beginning of
       * each timestep.
       */
      void update() override;

      /**
       * Evaluates values of material parameters and reaction terms.
       */
      void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                    typename Interface<dim>::MaterialModelOutputs &out) const override;

      /**
       * @name Reference quantities
       * @{
       */
      //double reference_viscosity() const override;

      double reference_darcy_coefficient() const override;

      /**
       * @}
       */

        /**
         * Compute the equilibrium melt fractions for the given input conditions.
         * @p in and @p melt_fractions need to have the same size.
         *
         * @param in Object that contains the current conditions.
         * @param melt_fractions Vector of doubles that is filled with the
         * equilibrium melt fraction for each given input conditions.
         * NOT IMPLEMENTED
         */
        void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                             std::vector<double> &melt_fractions) const override;

      /**
       * @name Functions used in dealing with run-time parameters
       * @{
       */
      /**
       * Declare the parameters this class takes through input files.
       */
      static void
      declare_parameters(ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter file.
       */
      void
      parse_parameters(ParameterHandler &prm) override;
      /**
       * @}
       */

      void
      create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const override;

    private:
      double reference_rho_s;
      double reference_rho_f;
      double reference_T;
      double eta_0;
      double epsdot_0;
      double eta_f;
      double activation_energy;
      double stress_exponent;
      double thermal_expansivity;
      double reference_specific_heat;
      double thermal_conductivity;
      double reference_permeability;
      double alpha_phi;
      double composition_density_change;
      bool include_melting_and_freezing;
      double melting_time_scale;
      double alpha_composition;
      double C_reference;
      double delta_eta_composition_max;

      double wBt;
      double wMu;
      double wMu0;
      double wMu1;
      double DT;
      
      double max_depth;

      /*
       * Calculate the equilibrium solid and melt composition
       */
      virtual void
      c_sf(const double temperature,
           const double pressure,
           double &c_s,
           double &c_f,
           int &PTfield=0) const;
    };

  }
}

#endif