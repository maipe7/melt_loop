/*
    Copyright (C) 2025 Petra Maierova

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

#include "melt_petrol.h"

#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>

#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>
#include <aspect/postprocess/melt_statistics.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    double
    MeltPetrol<dim>::
        reference_darcy_coefficient() const
    {
      double reference_porosity = 0.01; // 1% of melt
      return reference_permeability * std::pow(reference_porosity, 3.0) * std::pow(1.0 - reference_porosity, 2.0) / eta_f;
    }

    template <int dim>
    bool
    MeltPetrol<dim>::
        is_compressible() const
    {
      return false;
    }

    template <int dim>
    void
    MeltPetrol<dim>::
        melt_fractions(const MaterialModel::MaterialModelInputs<dim> &in,
                       std::vector<double> &melt_fractions) const
    // filled with zeroes - to be removed??
    {
      for (unsigned int q = 0; q < in.n_evaluation_points(); ++q)
      {
        melt_fractions[q] = 0.0;
      }
    }

    inline double fTwetsolid(double pressure)
    // pressure-dependent wet solidus temperature (~Holtz 2001) combined with ???;
    // works above 0.1 GPa
    // input pressure in Pa
    {
      const double pressure0 = 1.;
      pressure *= 1e-9;
      const double minwetsolid = 650. + 273.; 
      return (pressure < pressure0 ? 
        minwetsolid + 150. * std::pow(pressure - pressure0, 4.) / pressure0 : 
        minwetsolid +  50. * std::pow(pressure - pressure0, 2.) / pressure0);
    }

    inline double fTmus(double pressure)
    // muscovite dehydration line according to Thermocalc;
    // experiments put it higher (Patino Douce and Harris, 1998)
    // input pressure in Pa
    {
      return 620.0 + 273. + 130. * pressure * 1e-9;
    }

    template <int dim>
    void
    MeltPetrol<dim>::
        c_sf(const double temperature,
             const double pressure, double &c_s, double &c_f, int &PTfield) const
    {
      const double pressure0 = 1e9; // reference pressure for phase diagram construction

      // temperature of muscovite dehydration reaction:
      double Tmus = fTmus(pressure);
      // wet solidus temperature:
      double Twetsolid = fTwetsolid(pressure);

      // liquidus curve - fit to data in Makhluf et al., 2017 (based on Johannes, Holtz):
      // wliquidus(P,T)= (a*P+b)*(1-((T-c)/(d+e*P))**f) //P..GPa,T..C
      //const double c = 600. + 273., aG = 19e-9, b = 29., d = 360., eG = 205e-9, f = 0.13;
      const double c = 644. + 273., aG = 7.5*1e-9, b = 8.4, d = 316., eG = 198.0*1e-9, f = 0.41;
      double wliquidus = (aG * pressure + b) * (1 - std::pow((std::max(Twetsolid, temperature) - c) / (d + eG * pressure), f));
      // solidus temperature for c=0; dry solidus==dry liquidus - curves cross there:
      double Tdrysolid = c + d + eG * pressure;
      // liquid composition at temperature = wet solidus - DT; curves cross there:
      // TODO assert (Twetsolid - DT - c) > 0:
      AssertThrow(Twetsolid - DT - c > 0, ExcMessage("(Minimum wet solidus - Temperature width of dehydration reaction has to be > C in liquidus definition!"));
      double wmax = (aG * pressure + b) * (1 - std::pow((Twetsolid - DT - c) / (d + eG * pressure), f));

      // Construction of the solidus composition, a piecewise-linear function of temperature
      double T0 = Tdrysolid;
      double T1 = Tmus;
      double T2 = Tmus - DT;
      double T3 = Twetsolid;
      double T4 = Twetsolid - DT;
      double w0 = 0;
      double w1 = wBt;
      double w2 = wBt + wMu0;
      double w3 = wBt + wMu0 + wMu1;
      double w4 = wmax;
      double wlin1, wlin2, wlin3, wlin4;

      if (true) // pressure-dependent w1,2,3
      { 
        double w1Pvar = 0.0, w2Pvar=0.0, w3Pvar=0.0; // zero
        //double w1Pvar = -0.5, w2Pvar=0.0, w3Pvar=0.4; // pelitic // w1=0.5; w2=1.5; w3=2.0; w1Pvar=-0.5; w2Pvar=0.0; w3Pvar=0.4
        //double w1Pvar = -0.2, w2Pvar=0.0, w3Pvar=0.4; // granitic  // w1=0.2; w2=0.6; w3=1.0; w1Pvar=-0.2; w2Pvar=0.0; w3Pvar=0.4 # new
        w1 += w1Pvar * (pressure - pressure0) / pressure0;
        w2 += w2Pvar * (pressure - pressure0) / pressure0;
        w3 += w3Pvar * (pressure - pressure0) / pressure0;
      }
      wlin1 = w0 + (w1 - w0) * (temperature - T0) / (T1 - T0);
      wlin2 = w1 + (w2 - w1) * (temperature - T1) / (T2 - T1);
      wlin3 = w2 + (w3 - w2) * (temperature - T2) / (T3 - T2);
      wlin4 = w3 + (w4 - w3) * (temperature - T3) / (T4 - T3);
//cout << pressure << " " << temperature << " " << (Twetsolid - DT - c) << " " << (d + eG * pressure) << " " << w4 << " \n";
      double wsolidus =
          (temperature > Twetsolid - DT ? (temperature > Twetsolid ? (temperature > Tmus - DT ? (temperature > Tmus ? wlin1 : wlin2) : wlin3) : wlin4)
                                        : wliquidus);
      // fill in the output variables:
      if (temperature >= Tdrysolid) // above liquidus for all c
      {
        c_s = 0.0;
        c_f = 0.0;
        PTfield = 2;
      }
      else if (temperature < Twetsolid - DT) // below solidus for all c
      {
        c_s = wliquidus; // is the same as wsolidus
        c_f = wliquidus;
        PTfield = 0;
      }
      else
      {
        c_s = wsolidus;
        c_f = wliquidus;
        PTfield = 1;
      }
    }

    template <int dim>
    void
    MeltPetrol<dim>::
        evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      std::vector<double> old_porosity(in.n_evaluation_points());
      std::vector<double> old_peridotite(in.n_evaluation_points());
      std::vector<double> old_peridotiteF(in.n_evaluation_points());
      //std::vector<double> old_plastic_strain(in.n_evaluation_points());
      std::vector<double> old_total_strain(in.n_evaluation_points());
      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();
      if (this->include_melt_transport() && in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0 && !this->get_parameters().use_operator_splitting)
      {
        AssertThrow(false, ExcMessage("Melting model does not work without operator splitting."));
      }
      // Prepare reaction timestep for calculation of reaction terms:
      double dtt = 0;
      if (in.requests_property(MaterialProperties::reaction_terms))
      {
        const unsigned int number_of_reaction_steps = std::max(static_cast<unsigned int>(this->get_timestep() / this->get_parameters().reaction_time_step),
                                                               std::max(this->get_parameters().reaction_steps_per_advection_step, 1U));
        dtt = this->get_timestep() / static_cast<double>(number_of_reaction_steps);
      }
      /*
      const double reference_depth = 0.0;
      const Point<dim> representative_point = this->get_geometry_model().representative_point(reference_depth);
      const double reference_gravity = this->get_gravity_model().gravity_vector(representative_point).norm();
      */

      const bool is_rigid_present = this->introspection().compositional_name_exists("rigid");
      const bool is_totalstrain_present = this->introspection().compositional_name_exists("total_strain");

      for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
      { // TODO move this outside for loop?
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
        const unsigned int peridotiteF_idx = this->introspection().compositional_index_for_name("peridotiteF");
        //const unsigned int plastic_strain_idx = this->introspection().compositional_index_for_name("plastic_strain");
        int total_strain_idx = -1;
        if (is_totalstrain_present)
          total_strain_idx = this->introspection().compositional_index_for_name("total_strain");
        old_peridotite[i] = in.composition[i][peridotite_idx];
        old_peridotiteF[i] = in.composition[i][peridotiteF_idx]; //
        old_porosity[i] = in.composition[i][porosity_idx];
        //old_plastic_strain[i] = in.composition[i][plastic_strain_idx];
        if (is_totalstrain_present)
          old_total_strain[i] = in.composition[i][total_strain_idx];
        else
          old_total_strain[i] = 0;
          const double porosity = std::min(1.0, std::max(old_porosity[i], 0.0));
        // this choice (crop porosity or not) can make an important difference:
        double c_tot_old = old_peridotite[i] * (1.0 - porosity) + old_peridotiteF[i] * porosity;
        // double c_tot_old = old_peridotite[i] * (1.0 - old_porosity[i]) + old_peridotiteF[i] * old_porosity[i];
          
        // calculate strain-rate invariant:
        double edot_ii;
        const bool use_reference_strainrate = (this->get_timestep_number() == 0) &&
                                                (in.strain_rate[i].norm() <= std::numeric_limits<double>::min());
        if (use_reference_strainrate)
            edot_ii = epsdot_0;
        else
            edot_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(in.strain_rate[i])))), epsdot_0);
          
        // Calculate density:
        // temperature dependence of density is 1 - alpha * (T - Tref)
        double temperature_dependence = 1.0;
        temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
        // calculate composition dependence of density
        // TODO figure out correct way to treat it (based on observations/experiments):
        double delta_rho = this->introspection().compositional_name_exists("peridotite")
                                //   ? composition_density_change * (C_reference - c_tot_old) : 0.0;
                                 ? composition_density_change * (C_reference - old_peridotite[i]) : 0.0;
        delta_rho = 0.01 * std::max (std::min ( delta_rho, max_composition_density_change), -max_composition_density_change);
        out.densities[i] = reference_rho_s * (1.0 + delta_rho) * temperature_dependence;
        // Calculate viscosity:
        out.viscosities[i] = eta_0;
        double viscous_viscosity = out.viscosities[i];
        double vp_viscosity = viscous_viscosity;
        if (in.requests_property(MaterialProperties::viscosity)) // needed, because the strain_rate may not be filled otherwise
        {
          Assert(this->get_timestep_number() <= 1 || std::isfinite(in.strain_rate[i].norm()),
                 ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                            "not filled by the caller."));
          // viscosity reduction due to porosity:
          out.viscosities[i] *= std::max(exp(-alpha_phi * porosity), 1e-4); // ref
          // out.viscosities[i] *= exp(-alpha_phi * porosity); // exp
          /*  // coke // fits Costa/Keller well for 5 orders of magnitude decrease:
          double A = 0.05;
          double width = 0.26;
          out.viscosities[i] *= exp(-alpha_phi * porosity) * ((porosity > A) ? (porosity < A + width ? (std::exp(-width / (-porosity + A + width)) /
                                            (std::exp(-width / (-porosity + A + width)) + std::exp(-width / (width - (-porosity + A + width)))))
                                         : 0.0)
                                         : 1.0);
          */
          // normalized solid composition:
          const double C_solid_normalized = (C_reference - old_peridotite[i]);
          // const double C_solid_normalized = std::max(-1.0, std::min(1.0, (old_peridotite[i] - C_reference) / dC_solidus_liquidus));
          //  composition-dependent term:
          const double visc_composition_dependence = 
            std::max( std::min(exp(alpha_composition * C_solid_normalized), max_delta_eta_composition), 1./max_delta_eta_composition);
          out.viscosities[i] *= visc_composition_dependence;
          // temperature dependence of viscosity:
          double visc_temperature_dependence = 1.0;
          if (activation_energy != 0.0)
          {
            const double invtemperature = 1.0 / std::max(in.temperature[i], 1.0);
            visc_temperature_dependence = std::exp(activation_energy * invtemperature / stress_exponent / 8.31);
          }
          out.viscosities[i] *= visc_temperature_dependence;
          // stress (strain-rate) dependence of viscosity:
          out.viscosities[i] *= std::pow(edot_ii, ((1.0 - stress_exponent) / stress_exponent));
          // integrate strain-rate invariant to total strain:
          double e_ii = old_total_strain[i] + edot_ii*this->get_timestep();
          // strain-dependent weakening of ductile part:
          //out.viscosities[i] *= std::max(1.0-e_ii/e_ii_max*(1-visc_rel_drop),visc_rel_drop);
          out.viscosities[i] *= std::pow(10.0,std::max(-e_ii/e_ii_max*visc_rel_drop,-visc_rel_drop));
          viscous_viscosity = out.viscosities[i];
          vp_viscosity = viscous_viscosity;
          if(true && porosity < 1e-3 ) // Drucker-Prager TODO change
          {
            //double e_ii = old_plastic_strain[i] + edot_ii*this->get_timestep();
            double curr_friction_angle = friction_angle*std::max(1.0-e_ii/e_ii_max*(1-angle_rel_drop),angle_rel_drop);
            const double sin_phi = std::sin(curr_friction_angle);
            const double cos_phi = std::cos(curr_friction_angle);
            const double yield_stress = std::min(cohesion * cos_phi + std::max(0.0,in.pressure[i]) * sin_phi,1e12);
            double current_stress = 2. * viscous_viscosity * edot_ii;
            if (current_stress >= yield_stress) {
              vp_viscosity = yield_stress*0.5/edot_ii;
              //std::cout << " " << yield_stress << " " << current_stress << " " << edot_ii << " " << out.viscosities[i] << "\n";
            }
            // eta_eff = min ( eta_eff, yield/2.0/edot_ii ); the same ...
          }
          // prescribe maximum viscosity where "rigid" field is present:
          if (is_rigid_present)
          {
            unsigned int rigid_idx = this->introspection().compositional_index_for_name("rigid");
            double rigid = std::min(std::max(in.composition[i][rigid_idx], 0.0), 1.0);
            vp_viscosity = vp_viscosity + (1e23-vp_viscosity)*rigid;
          } 
          // total value of shear viscosity is cropped to min and max values
          out.viscosities[i] = std::min(std::max(vp_viscosity, 1e16), 1e23); // TODO modif 1e16

        }
        // Fill in melt material properties:
        if (melt_out != nullptr)
        {
          melt_out->permeabilities[i] = reference_permeability * std::pow(porosity, 3) * std::pow(1.0 - porosity, 2);
          melt_out->fluid_density_gradients[i] = Tensor<1, dim>(); // not calculated

          melt_out->fluid_viscosities[i] = eta_f;
          if (!use_constant_eta_f)
          {
            double waterWtPc = std::max(0.1, old_peridotiteF[i]);
            if (false) //
            {
              /*function to fit temperature dependence of data (Holz 01, Fig. 4):
              A=174.; B=-0.015; T0=16.6; P0=8.; dWdP=0.25; Tinfty=500.
              Wfluid(T,P)=A/(T-Tinfty)+B*T+T0+(P-P0)*dWdP*/
              const double temperatureReduced = std::max(in.temperature[i], 823.);
              waterWtPc = 174. / (temperatureReduced - 773.) - 0.015 * (temperatureReduced - 273.) + 16.6 + (in.pressure[i] / 1e8 - 8.) * 0.25;
            }
            const double invtemperature = 1.0 / std::max(in.temperature[i], 1.0);
            /* viscosity according to Schulze etal 1995:
             LogEta0=-2.5726 #-1.5726 for equation in Poise
             EA(W)=(448.03-252.12*W**0.11)*1e3
             LogEta(T,W)=LogEta0+EA(W)/(2.303*R*T)
             eta(T,W)=10.0**LogEta(T,W) */
            double eta_f;
            if (melt_viscosity_law==2)
              eta_f = std::pow(10.0, -7.5461 + 16280.*invtemperature + (0.59784-1235.4*invtemperature)*waterWtPc); // too low?! // Scaillet
            else
              eta_f = std::pow(10.0, -2.5726 + (448.03 - 252.12 * std::pow(waterWtPc, 0.11)) * 1e3 / (2.303 * 8.31) * invtemperature); // Schulze
            melt_out->fluid_viscosities[i] = std::max(std::min(eta_f, 1e7), 1e2); // TODO test without cropping
          }

          // Calculate melt density:
          // TODO find correct way:
          //const double delta_rho = 0.0; compositional dependence of melt density is neglected
          delta_rho = this->introspection().compositional_name_exists("peridotiteF")
                                ? composition_density_change * (C_reference - old_peridotiteF[i]) : 0.0;
          delta_rho = 0.01 * std::max (std::min ( delta_rho, max_composition_density_change), -max_composition_density_change);

          double temperature_dependence = 1.0;
          temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
          melt_out->fluid_densities[i] = reference_rho_f * (1.0 + delta_rho) * temperature_dependence;
          // compaction viscosity defined relatively to shear viscosity
          melt_out->compaction_viscosities[i] = out.viscosities[i] * (2.0 + 0.3 / std::max(porosity, 3e-4)); // ref
          // melt_out->compaction_viscosities[i] = out.viscosities[i] / (1.0 - std::min(porosity, 1.0-3e-4))/ (std::max(porosity, 3e-4)); // 1-phi
          // melt_out->compaction_viscosities[i] = out.viscosities[i] / (std::max(porosity, 3e-4)); // 1over
          // melt_out->compaction_viscosities[i] = out.viscosities[i]*(2.0-2.7*std::log10(std::max(porosity,1e-3))); // uncomment for logarithmic dependence of bulk viscosity
        }

        // Calculate melting/freezing and composition reactions:
        // By default, no melting or freezing --> set all reactions to zero
        for (unsigned int c = 0; c < in.composition[i].size(); ++c)
        {
          out.reaction_terms[i][c] = 0.0;
          if (this->get_parameters().use_operator_splitting && reaction_rate_out != nullptr)
            reaction_rate_out->reaction_rates[i][c] = 0.0;
        }

        if (include_melting_and_freezing && in.requests_property(MaterialProperties::reaction_terms))
        { 
          //  Calculate solid and fluid compositions:
          double c_s;
          double c_f;
          int PTfield = 0;
          //double lithostatic_pressure = reference_gravity * reference_rho_s * (max_depth - in.position[i][1]); //  (possible modification with depth = this->get_geometry_model().depth(point);) // works only without mesh deformation
          //c_sf(in.temperature[i], lithostatic_pressure, c_s, c_f, PTfield); // use lithostatic pressure - mesh_deformation approximated by max topography (see function update())
          double dynamic_pressure = std::max(1.0,in.pressure[i]); // pressure on top set in BCs, but always positive
          c_sf(in.temperature[i], dynamic_pressure, c_s, c_f, PTfield); // use dynamic pressure

          if (c_f < c_s + 1e-5)
            c_f = c_s + 1e-5; // avoid possible division by zero and switch between c_f and c_s

          bool is_rigid = false;
          if (is_rigid_present)
          {
            unsigned int rigid_idx = this->introspection().compositional_index_for_name("rigid");
            is_rigid = in.composition[i][rigid_idx] > 0.01;
          } 

          // Calculate updates of porosity and compositions.
          // Updates of solid and melt compositions are calculated from mass conservation
          // assuming equal densities of solid and melt.
          for (unsigned int c = 0; c < in.composition[i].size(); ++c)
          { // fill reaction rate outputs as the model uses operator splitting
            if (this->get_parameters().use_operator_splitting && this->get_timestep_number() > 0)
            {
              if (reaction_rate_out != nullptr)
              {
                if (c_tot_old < c_s || is_rigid) //... below solidus - equilibrium porosity is 0
                {
                  if (c == porosity_idx)
                    reaction_rate_out->reaction_rates[i][c] =
                        -old_porosity[i] / melting_time_scale;
                  else if (c == peridotite_idx || c == peridotiteF_idx)
                    reaction_rate_out->reaction_rates[i][c] =
                          old_porosity[i] * (old_peridotiteF[i] - old_peridotite[i]) / melting_time_scale;
                 /* else if (c == peridotite_idx)
                    if (PTfield > 0) // || porosity > 1e-3) // TODO rethink usage of PTfield, TODO epsilon
                      reaction_rate_out->reaction_rates[i][c] =
                          porosity * (old_peridotiteF[i] - old_peridotite[i]) / melting_time_scale;
                    else // below wet solidus - solid composition tends to bulk composition
                      reaction_rate_out->reaction_rates[i][c] =
                          (c_tot_old - old_peridotite[i]) / melting_time_scale;
                  else if (c == peridotiteF_idx)
                    if (PTfield > 0) //|| porosity > 1e-3) // TODO rethink usage of PTfield, TODO epsilon
                      reaction_rate_out->reaction_rates[i][c] =
                          // porosity * (old_peridotiteF[i] - old_peridotite[i]) / melting_time_scale; //
                          (c_f - old_peridotiteF[i]) / melting_time_scale; // adjust to liquidus composition // TODO which one is correct? Why??
                    else                                                   // below wet solidus - liquid composition arbitrary, artificially adjusted to liquidus composition
                      reaction_rate_out->reaction_rates[i][c] =
                          (c_f - old_peridotiteF[i]) / melting_time_scale; */
                  else
                    reaction_rate_out->reaction_rates[i][c] = 0.0;
                }

                else if (c_tot_old < c_f) // ... between solidus and liquidus
                { //std::cout << c_tot_old << " " << c_f << "above solidus \n";
                  const double Dc_s = c_s - old_peridotite[i];
                  const double Dc_f = c_f - old_peridotiteF[i];
                  // porosity update derived from lever rule
                  const double Dpor = ((1.0 - old_porosity[i]) * Dc_s + old_porosity[i] * Dc_f) /
                                      ((old_peridotite[i] + Dc_s / melting_time_scale * dtt) - (old_peridotiteF[i] + Dc_f / melting_time_scale * dtt));
                  if (c == porosity_idx)
                    reaction_rate_out->reaction_rates[i][c] = Dpor / melting_time_scale;
                  else if (c == peridotite_idx)
                    reaction_rate_out->reaction_rates[i][c] =
                        Dc_s / melting_time_scale;
                  else if (c == peridotiteF_idx)
                    reaction_rate_out->reaction_rates[i][c] =
                        Dc_f / melting_time_scale;
                  else
                    reaction_rate_out->reaction_rates[i][c] = 0.0;
                }
                else // if ( c_tot_old > c_f) // temperature (composition) above liquidus - equilibrium porosity is 1
                {    // TODO porosity-old_porosity?
                   //std::cout << c_tot_old << " " << c_f <<  "above liquidus \n";
                  if (c == porosity_idx)
                    reaction_rate_out->reaction_rates[i][c] =
                        (1.0 - old_porosity[i]) / melting_time_scale;
                  else if (c == peridotite_idx || c == peridotiteF_idx)
                    reaction_rate_out->reaction_rates[i][c] =
                        (1.0 - old_porosity[i]) * (old_peridotite[i] - old_peridotiteF[i]) / melting_time_scale;
                  /*else if (c == peridotite_idx)
                    reaction_rate_out->reaction_rates[i][c] =
                        //(c_tot_old - old_peridotite[i]) / melting_time_scale; // ??
                        (1.0 - porosity) * (old_peridotite[i] - old_peridotiteF[i]) / melting_time_scale; // may produce negative cs if we start out of equilibrium
                  else if (c == peridotiteF_idx)
                    reaction_rate_out->reaction_rates[i][c] =
                        //(1.0 - porosity) * (old_peridotite[i] - old_peridotiteF[i]) / melting_time_scale;
                        (c_tot_old - old_peridotiteF[i]) / melting_time_scale; // adjust liquid composition to bulk composition
                        */
                  else
                    reaction_rate_out->reaction_rates[i][c] = 0.0;
                }
                //if (c == plastic_strain_idx)
                //  reaction_rate_out->reaction_rates[i][c] = edot_ii*(1.0-vp_viscosity/viscous_viscosity); // TODO check division by zero
              }
              out.reaction_terms[i][c] = 0.0;
              if (c == total_strain_idx) 
              {
                out.reaction_terms[i][c] = this->get_timestep()*edot_ii; // TODO better check division by zero
              }
              /*if (c == plastic_strain_idx) 
              {
                // if porosity>0 relax plastic strain to 0
                if (old_porosity[i]>0.01)
                  out.reaction_terms[i][c] = - old_plastic_strain[i]*this->get_timestep()*3e-14*old_porosity[i];  // TODO improve to be time-step independent
                else
                  out.reaction_terms[i][c] = this->get_timestep()*edot_ii*(1.0-vp_viscosity/std::max(1.0,viscous_viscosity)); // TODO better check division by zero
              }*/
            }
          }
        }

        // a hack to account for melt buoyancy even in models without melt migration:
        if (!this->include_melt_transport()) {
          out.densities[i] *= 1.0 - porosity*(1.0 - reference_rho_f/reference_rho_s);
        } 
        out.entropy_derivative_pressure[i] = 0.0;
        out.entropy_derivative_temperature[i] = 0.0;
        out.thermal_expansion_coefficients[i] = thermal_expansivity;
        out.specific_heat[i] = reference_specific_heat;
        out.thermal_conductivities[i] = thermal_conductivity;
        out.compressibilities[i] = 0.0;
      }
    }

    template <int dim>
    void
    MeltPetrol<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt petrol");
        {
          prm.declare_entry("Reference solid density", "3000.",
                            Patterns::Double(0.),
                            "Reference density of the solid $\\rho_{s,0}$. "
                            "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry("Reference melt density", "2500.",
                            Patterns::Double(0.),
                            "Reference density of the melt$\\rho_{f,0}$. "
                            "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry("Reference temperature", "293.",
                            Patterns::Double(0.),
                            "The reference temperature $T_0$. The reference temperature is used "
                            "in the density calculation. Units: \\si{\\kelvin}.");
          prm.declare_entry("Reference shear viscosity", "5e20",
                            Patterns::Double(0.),
                            "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                            "This viscosity may be modified by temperature, strain-rate, composition and porosity "
                            "dependencies. Units: \\si{\\pascal\\second}.");
          prm.declare_entry("Reference strainrate", "1e-14",
                            Patterns::Double(0.),
                            "Reference value of strain rate for viscosity calculation. Used "
                            "as a minimum value of strain rate."
                            "Units: \\si{\\per\\second}.");
          prm.declare_entry("Stress viscosity exponent", "1.",
                            Patterns::Double(0.),
                            "Stress exponent for viscosity calculation.");
          prm.declare_entry("Reference melt viscosity", "0.",
                            Patterns::Double(),
                            "The value of the reference melt viscosity $\\eta_f$ "
                            "used for reference Darcy coefficient calculation."
                            "Units: \\si{\\pascal\\second}.");
          prm.declare_entry("Use constant melt viscosity", "true",
                            Patterns::Bool(),
                            "Whether to use constant melt viscosity or rather "
                            "composition (interpreted as wt% water-content)- dependent one.");
          prm.declare_entry("Exponential melt weakening factor", "27.",
                            Patterns::Double(0.),
                            "The porosity dependence of the solid viscosity. Units: dimensionless.");
          prm.declare_entry("Activation energy", "0.0",
                            Patterns::Double(0.),
                            "Activation energy of the solid viscosity. Units: \\si{\\joule\\per\\mole}");
          prm.declare_entry("Thermal conductivity", "4.7",
                            Patterns::Double(0.),
                            "The value of the thermal conductivity $k$. "
                            "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry("Reference specific heat", "1250.",
                            Patterns::Double(0.),
                            "The value of the specific heat $C_p$. "
                            "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry("Thermal expansion coefficient", "2e-5",
                            Patterns::Double(0.),
                            "The value of the thermal expansion coefficient $\\beta$. "
                            "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry("Reference permeability", "1e-8",
                            Patterns::Double(),
                            "Reference permeability of the solid host rock."
                            "Units: \\si{\\meter\\squared}.");
          prm.declare_entry("Composition pc density change", "0.0",
                            Patterns::Double(),
                            "Change of density due to compositional change of 1 (in percent). "
                            "Units: percent.");
          prm.declare_entry("Maximum composition pc density change", "0.0",
                            Patterns::Double(),
                            "Maximum change of density due to composition. "
                            "Units: percent.");
          prm.declare_entry("Surface solidus", "1300.",
                            Patterns::Double(0.),
                            "Solidus at the surface (zero lithostatic pressure). "
                            "Units: \\si{\\kelvin}.");
          prm.declare_entry("Solidus liquidus difference", "500.",
                            Patterns::Double(),
                            "Difference between solidus and liquidus temperatures. "
                            "Units: \\si{\\kelvin}.");
          prm.declare_entry("Solidus liquidus composition difference", "1.",
                            Patterns::Double(),
                            "Difference between solidus and liquidus compositions. "
                            "Units: ---.");
          prm.declare_entry("Reference solid composition", "0.",
                            Patterns::Double(0.),
                            "Reference composition of solid. "
                            "Units: ---.");
          prm.declare_entry("Pressure solidus change", "0.0",
                            Patterns::Double(),
                            "The linear solidus temperature change with pressure. For positive "
                            "values, the solidus gets increased for positive pressures. "
                            "Units: \\si{\\kelvin\\per\\pascal}.");
          prm.declare_entry("Include melting and freezing", "true",
                            Patterns::Bool(),
                            "Whether to include melting and freezing (according to a simplified "
                            "linear melting approximation in the model (if true), or not (if "
                            "false).");
          prm.declare_entry("Melting time scale for operator splitting", "1e3",
                            Patterns::Double(0.),
                            "In case the operator splitting scheme is used, the porosity field can not "
                            "be set to a new equilibrium melt fraction instantly, but the model has to "
                            "provide a melting time scale instead. This time scale defines how fast melting "
                            "happens, or more specifically, the parameter defines the time after which "
                            "the deviation of the porosity from the equilibrium melt fraction will be "
                            "reduced to a fraction of $1/e$. So if the melting time scale is small compared "
                            "to the time step size, the reaction will be so fast that the porosity is very "
                            "close to the equilibrium melt fraction after reactions are computed. Conversely, "
                            "if the melting time scale is large compared to the time step size, almost no "
                            "melting and freezing will occur."
                            "\n\n"
                            "Also note that the melting time scale has to be larger than or equal to the reaction "
                            "time step used in the operator splitting scheme, otherwise reactions can not be "
                            "computed. If the model does not use operator splitting, this parameter is not used. "
                            "Units: yr or s, depending on the ``Use years "
                            "in output instead of seconds'' parameter.");
          prm.declare_entry("Exponential compositional viscosity factor", "0.0",
                            Patterns::Double(0.),
                            "Exponential dependency of solid viscosity on the "
                            "solid composition (peridotite field). "
                            "Dimensionless factor. With a value of 0.0 (the default) the "
                            "viscosity does not depend on the composition.");
          prm.declare_entry("Maximum rel compositional viscosity change", "1.0",
                            Patterns::Double(0.),
                            "Maximum relative change of solid viscosity due to composition variations.");

          prm.declare_entry("Water in biotite", "0.5",
                            Patterns::Double(0.),
                            "H2O in rock bond in biotite (Wt%).");
          prm.declare_entry("Water in muscovite", "0.5",
                            Patterns::Double(0.),
                            "Total H2O in rock bond in muscovite (Wt%).");
          prm.declare_entry("Jump in muscovite", "0.2",
                            Patterns::Double(0.),
                            "How much water is released during muscovite dehydration reaction (in Wt%).");
          prm.declare_entry("Temperature width of reactions", "5.",
                            Patterns::Double(0.),
                            "Width of the temperature interval, over which the abrupt dehydration reaction and wet solidus are smoothed (K).");
          prm.declare_entry("Melt viscosity law", "1",
                            Patterns::Integer(0),
                            "Define law for melt viscosity:"
                            "1 = Schulze et al. 1996, default; 2 = Scaillet et al. 1996.");
          prm.declare_entry("Angle of internal friction", "10.",
                            Patterns::Double(0.),
                            "Angle of internal friction (degrees).");
                            prm.declare_entry("Cohesion", "1e30",
                              Patterns::Double(0.),
                              "Cohesion (Pa).");
          prm.declare_entry("Maximum strain for weakening", "1.0",
                            Patterns::Double(0.),
                              "Maximum strain for weakening of friction angle.");
          prm.declare_entry("Friction angle relative weakening", "1.0",
                            Patterns::Double(0.),
                              "Final relative weakening of friction angle. Default no weakening.");
          /*prm.declare_entry("Relative ductile weakening", "1.0",
                            Patterns::Double(0.),
                               "Final relative weakening of ductile part of viscosity. Default no weakening.");*/
          prm.declare_entry("Relative log10 ductile weakening", "0.0",
                                Patterns::Double(0.),
                                   "Final relative weakening of ductile part of viscosity in logscale. Default no weakening.");
    
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MeltPetrol<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt petrol");
        {
          reference_rho_s = prm.get_double("Reference solid density");
          reference_rho_f = prm.get_double("Reference melt density");
          reference_T = prm.get_double("Reference temperature");
          eta_0 = prm.get_double("Reference shear viscosity");
          epsdot_0 = prm.get_double("Reference strainrate");
          eta_f = prm.get_double("Reference melt viscosity");
          use_constant_eta_f = prm.get_bool("Use constant melt viscosity");
          reference_permeability = prm.get_double("Reference permeability");
          activation_energy = prm.get_double("Activation energy");
          stress_exponent = prm.get_double("Stress viscosity exponent");
          friction_angle = prm.get_double("Angle of internal friction")*constants::degree_to_radians; // degrees -> rad;
          cohesion = prm.get_double("Cohesion");
          e_ii_max = prm.get_double("Maximum strain for weakening"); // TODO assert positivity
          angle_rel_drop = prm.get_double("Friction angle relative weakening"); // TODO assert positivity
          visc_rel_drop = prm.get_double("Relative log10 ductile weakening"); // TODO assert positivity
          thermal_conductivity = prm.get_double("Thermal conductivity");
          reference_specific_heat = prm.get_double("Reference specific heat");
          thermal_expansivity = prm.get_double("Thermal expansion coefficient");
          alpha_phi = prm.get_double("Exponential melt weakening factor");
          composition_density_change = prm.get_double("Composition pc density change");
          max_composition_density_change = prm.get_double("Maximum composition pc density change");

          wBt = prm.get_double("Water in biotite");
          wMu = prm.get_double("Water in muscovite");
          wMu0 = prm.get_double("Jump in muscovite");
          wMu1 = wMu - wMu0;
          DT = prm.get_double("Temperature width of reactions");

          C_reference = prm.get_double("Reference solid composition");
          include_melting_and_freezing = prm.get_bool("Include melting and freezing");
          melting_time_scale = prm.get_double("Melting time scale for operator splitting");
          alpha_composition = prm.get_double("Exponential compositional viscosity factor");
          max_delta_eta_composition = prm.get_double("Maximum rel compositional viscosity change");
          melt_viscosity_law = prm.get_integer("Melt viscosity law");

          if (melt_viscosity_law!=1 && melt_viscosity_law!=2)
            AssertThrow(false, ExcMessage("Unknown melt viscosity law!"));

          AssertThrow(stress_exponent > 0, ExcMessage("Stress exponent has to be > 0!"));

          AssertThrow(epsdot_0 > 0, ExcMessage("Reference strainrate has to be > 0!"));

          if (this->convert_output_to_years() == true)
            melting_time_scale *= year_in_seconds;

          if (this->include_adiabatic_heating() == true)
          {
            AssertThrow(false,
                        ExcMessage("Adiabatic heating not implemented."));
          }

          // assert operator splitting
          if (this->get_parameters().use_operator_splitting)
          {
            AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                        ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step) + " in the operator splitting scheme is too large to compute melting rates! "
                                                                                                                                 "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                                                                                                                 "operator splitting' chosen in the material model, which is currently " +
                                   Utilities::to_string(melting_time_scale) + "."));
            AssertThrow(melting_time_scale > 0,
                        ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));
            AssertThrow(this->introspection().compositional_name_exists("porosity"),
                        ExcMessage("Material model Melt petrol with melt transport only "
                                   "works if there is a compositional field called porosity."));
          }

          // if (this->include_melt_transport())
          {
            AssertThrow(this->get_parameters().use_operator_splitting,
                        ExcMessage("Material model Melt petrol with melt transport only "
                                   "works with operator splitting."));
          }
          AssertThrow(this->introspection().compositional_name_exists("porosity"),
                      ExcMessage("Material model Melt petrol only works if there is a "
                                 "compositional field called porosity."));
          AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                      ExcMessage("Material model Melt petrol only works if there is a "
                                 "compositional field called peridotite."));
          AssertThrow(this->introspection().compositional_name_exists("peridotiteF"),
                      ExcMessage("Material model Melt petrol only works if there is a "
                                 "compositional field called peridotiteF."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MeltPetrol<dim>::create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().use_operator_splitting && out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr)
      {
        const unsigned int n_points = out.n_evaluation_points();
        out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>>(n_points, this->n_compositional_fields()));
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltPetrol,
                                   "melt petrol",
                                   "A material model for the modeling of melt transport. "
                                   "It includes a simple melting parametrization with "
                                   "piecewise-linear solidus and experiment-based liquidus temperatures, "
                                   "and equilibrium melt fraction calculated by the lever rule.")
  }
}
