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

#include "melt_global_mod.h"

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
        /*template <int dim>
        double
        MeltGlobalMod<dim>::
            reference_viscosity() const
        {
          return eta_0;
        }*/

    template <int dim>
    double
    MeltGlobalMod<dim>::
        reference_darcy_coefficient() const
    {
      double reference_porosity = 0.01; // 1% of melt
      return reference_permeability * std::pow(reference_porosity, 3.0) * std::pow(1.0 - reference_porosity, 2.0) / eta_f;
    }

    template <int dim>
    bool
    MeltGlobalMod<dim>::
        is_compressible() const
    {
      return false;
    }

    template <int dim>
    void
    MeltGlobalMod<dim>::
        melt_fractions(const MaterialModel::MaterialModelInputs<dim> &in,
                       std::vector<double> &melt_fractions) const
    {
      for (unsigned int q = 0; q < in.n_evaluation_points(); ++q)
      {
        melt_fractions[q] = 0.0;
      }
    }
/*

    template <int dim>
    double
    MeltGlobalMod<dim>::
    c_fluid (const double temperature,
             const double pressure) const
    {
      double T0=T00; //+ pressure_solidus_change*pressure;
      double Tlim=Tlim0; //+ pressure_solidus_change*pressure;
      double T05=1023.0; //or T0*0.25+Tlim*0.75;
      double CWmaxP=10.0+(pressure-1.0e9)*8e-9;

      double al=(-4.0*T05+2.0*(T0+Tlim))/(CWmaxP*CWmaxP);
      //TODO Assert(al<>0)
      double bl=(4.0*T05-Tlim-3.0*T0)/CWmaxP;
      double cl=T0;
      double disc= bl*bl-4.0*al*(cl-temperature);
      if ( temperature >= T0)
      {
        return 0.0;
      }
      else if (disc>0)
      {
        double C = (-bl - std::sqrt(disc) )/(2.0*al);
        //return 1.0 - C/CWmaxP;
        return C; // return absolute composition
      }
      else
        return CWmaxP; // TODO cbulk
    }


    template <int dim>
    double
    MeltGlobalMod<dim>::
    c_solid (const double temperature,
             const double pressure) const
    {
      double CWmaxP=10.0+(pressure-1.0e9)*8e-9;

      double T0=T00; // + pressure_solidus_change*pressure;
      double Tlim=Tlim0; // + pressure_solidus_change*pressure;
      double cs=Tlim-Tknee;
      double as=CWmaxP/(T0-Tlim)*(T0-cs)*(Tlim-cs); // TODO assert T0<>Tlim
      double bs=as/(T0-cs); // TODO assert T0 <> cs

      if (temperature >= T0 )
      {
        return 0.0;
      }
      else if (temperature >= Tlim )
      {
        double C = as/(temperature-cs)-bs;
        //return 1.0 - C/CWmaxP; // return scaled composition C'
        return C; // return absolute composition
      }
      else
        return CWmaxP; // TODO cbulk??  // TODO assert T<>cs??
    }
*/
/*
    template <int dim>
    double
    MeltGlobalMod<dim>::
    c_solid (const double temperature,
             const double pressure) const
    {
      double Tmax=960.+205*pressure*1e-9+273.;
      double Tlim=640.0+273.;
      double DTknee=10.0;
      double Tasymp=Tlim-DTknee;
      double CWmaxP=30.;
      double as=CWmaxP/(Tmax-Tlim)*(Tmax-Tasymp)*DTknee;
      double bs=as/(Tmax-Tasymp);
      if (temperature >= Tmax ) // above liquidus
      {
        return 0.0;
      }
      else if (temperature >= Tlim ) // above solidus
      {
        return as/(temperature-Tasymp)-bs;
      }
      else // below solidus
        return CWmaxP; // TODO cbulk??
    }

    template <int dim>
    double
    MeltGlobalMod<dim>::
    c_fluid (const double temperature,
             const double pressure) const
    {
      double c=600.+273.;
      double a=19.;
      double b=29.;
      double d=360.;
      double e=205.;
      double f=0.13;
      double Tmax= c+d+e*pressure*1e-9;

      double Tlim=640.0+273.;
      double DTknee=10.0;
      double Tasymp=Tlim-DTknee;
      double CWmaxP=30.;
      double as=CWmaxP/(Tmax-Tlim)*(Tmax-Tasymp)*DTknee;
      double bs=as/(Tmax-Tasymp);
      double csolidus=as/(temperature-Tasymp)-bs; // TODO div 0
      double cliquidus=(a*pressure*1e-9+b)*(1-std::pow((temperature-c)/(d+e*pressure*1e-9),f));
      
      //w(P,T)= (a*P+b)*(1-((T-c)/(d+e*P))**f) //P..GPa,T..C
      if ( temperature >= Tmax) // above liquidus
      {
        return 0.0;
      }
      else if (csolidus < cliquidus)
      {
        return cliquidus; // return absolute composition
      }
      else
        return CWmaxP; // TODO cbulk
    }
*/
    template <int dim>
    void
    MeltGlobalMod<dim>::
    c_sf (const double temperature,
          const double pressure, double &c_s, double &c_f) const
    {
      double c=600.+273.;
      //double a=19.;
      double aG=19e-9;
      double b=29.;
      double d=360.;
      //double e=205.;
      double eG=205e-9;
      double f=0.13;


      double Tmax=c+d+eG*pressure;
      //double Tmax= c+d+e*1e-9*pressure;

      //double Tlim=600.0+ 100.*std::pow((1e-9*pressure-1.),2) + 273.; // quadratic
      double Tlim;
      if (pressure > 1e9) Tlim=600.0+ 75.*(1e-9*pressure-1.) + 273.; // linear
      else Tlim=600.0+273.;
      double DTknee=2.0;
      double CWmax=30.;

      double Tasymp=Tlim-DTknee;
      double as=CWmax/(Tmax-Tlim)*(Tmax-Tasymp)*DTknee;
      double bs=as/(Tmax-Tasymp);
      double csolidus;
      double cliquidus;
      if (temperature>Tasymp) csolidus=as/(temperature-Tasymp)-bs;
      else csolidus=CWmax;
      cliquidus=std::min((aG*pressure+b)*(1-std::pow((temperature-c)/(d+eG*pressure),f)),CWmax);
      //cliquidus=std::min((1e-9*a*pressure+b)*(1-std::pow((temperature-c)/(d+1e-9*e*pressure),f)),CWmax);
      
      //w(P,T)= (a*P+b)*(1-((T-c)/(d+e*P))**f) //P..GPa,T..C
      if ( temperature >= Tmax) // above liquidus
      {
        c_s=0.0;
        c_f=0.0;
      }
      else if (csolidus < cliquidus)
      {
        c_s=csolidus;
        c_f=cliquidus;
      }
      else
      { 
        c_s=CWmax;
        c_f=CWmax; // TODO cbulk
      }
    }

    template <int dim>
    void
    MeltGlobalMod<dim>::
        evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      std::vector<double> old_porosity(in.n_evaluation_points());
      std::vector<double> old_peridotite(in.n_evaluation_points());
      std::vector<double> old_peridotiteF(in.n_evaluation_points());
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
      // Prepare gravity for calculation of lithostatic pressure:
      const double max_depth = this->get_geometry_model().maximal_depth();
      const double reference_depth = 0.0;
      const Point<dim> representative_point = this->get_geometry_model().representative_point(reference_depth);
      const double reference_gravity = this->get_gravity_model().gravity_vector(representative_point).norm();

      for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
        const unsigned int peridotiteF_idx = this->introspection().compositional_index_for_name("peridotiteF");
        old_peridotite[i] = in.composition[i][peridotite_idx];
        old_peridotiteF[i] = in.composition[i][peridotiteF_idx]; // 
        old_porosity[i] = in.composition[i][porosity_idx];
        const double porosity = std::min(1.0, std::max(old_porosity[i], 0.0));

        // Calculate density:
        // temperature dependence of density is 1 - alpha * (T - Tref)
        double temperature_dependence = 1.0;
        temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
        // calculate composition dependence of density
        const double delta_rho = this->introspection().compositional_name_exists("peridotite")
                                     ? composition_density_change * std::max(-1.0, std::min(1.0, (old_peridotite[i] - C_reference) / dC_solidus_liquidus))
                                     : 0.0; // delta_rho is cropped - peridotite only in <-1,1> is used.
        out.densities[i] = (reference_rho_s + delta_rho) * temperature_dependence;
        // Calculate viscosity:
        out.viscosities[i] = eta_0;
        if (in.requests_property(MaterialProperties::viscosity)) // needed, because the strain_rate may not be filled otherwise
        {
          Assert(this->get_timestep_number() <= 1 || std::isfinite(in.strain_rate[i].norm()),
                 ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                            "not filled by the caller."));
          // viscosity reduction due to porosity:
          out.viscosities[i] *= std::max(exp(-alpha_phi * porosity), 1e-4);
          // normalized solid composition:
          const double C_solid_normalized = std::max(-1.0, std::min(1.0, (old_peridotite[i] - C_reference) / dC_solidus_liquidus));
          // composition-dependent term:
          const double visc_composition_dependence = std::min(exp(alpha_composition * C_solid_normalized), delta_eta_composition_max);
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
          const bool use_reference_strainrate = (this->get_timestep_number() == 0) &&
                                                (in.strain_rate[i].norm() <= std::numeric_limits<double>::min());
          double edot_ii;
          if (use_reference_strainrate)
            edot_ii = epsdot_0;
          else
            edot_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(in.strain_rate[i])))), epsdot_0);
          out.viscosities[i] *= std::pow(edot_ii, ((1.0 - stress_exponent) / stress_exponent));
          // total value of shear viscosity is cropped to min and max values
          out.viscosities[i] = std::min(std::max(out.viscosities[i], 1e16), 1e23);
        }

        // Fill in melt material properties:
        if (melt_out != nullptr)
        {
          melt_out->fluid_viscosities[i] = eta_f;
          melt_out->permeabilities[i] = reference_permeability * std::pow(porosity, 3) * std::pow(1.0 - porosity, 2);
          melt_out->fluid_density_gradients[i] = Tensor<1, dim>(); // not calculated

          // Calculate melt density:
          // compositional dependence of melt density is neglected, as it is smaller than for solid density
          const double delta_rho = 0.0;
          double temperature_dependence = 1.0;
          temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
          melt_out->fluid_densities[i] = (reference_rho_f + delta_rho) * temperature_dependence;
          // compaction viscosity defined relatively to shear viscosity
          melt_out->compaction_viscosities[i] = out.viscosities[i] * (2.0 + 0.3 / std::max(porosity, 3e-4));
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

        //if (this->include_melt_transport())
        {
          if (include_melting_and_freezing && in.requests_property(MaterialProperties::reaction_terms))
          {
            double c_tot_old = old_peridotite[i] * (1.0 - old_porosity[i]) + old_peridotiteF[i] * old_porosity[i];
            // Calculate solid and fluid compositions:
            // (possible modification with depth = this->get_geometry_model().depth(point);)
            double lithostatic_pressure = reference_gravity * reference_rho_s * (max_depth - in.position[i][1]);
            //double c_s = c_solid(in.temperature[i], lithostatic_pressure);
            //double c_f = c_fluid(in.temperature[i], lithostatic_pressure);
            double c_s;
            double c_f;
            c_sf(in.temperature[i], lithostatic_pressure, c_s, c_f);
            
            if (c_f < c_s + 1e-5)
              c_f = c_s + 1e-5; // avoid (im)possible division by zero and switch between c_f and c_s

            // Calculate updates of porosity and compositions.
            // Updates of solid and melt compositions are calculated from mass conservation
            // assuming equal densities of solid and melt.
            for (unsigned int c = 0; c < in.composition[i].size(); ++c)
            { // fill reaction rate outputs as the model uses operator splitting
              if (this->get_parameters().use_operator_splitting && this->get_timestep_number() > 0)
              {
                if (reaction_rate_out != nullptr)
                {
                  //if (in.position[i][1]<10e3) cout << c_tot_old << " " << c_s << " " << c_f << " " << lithostatic_pressure << " " << in.temperature[i] ;

                  if (c_tot_old < c_s) //... below solidus - equilibrium porosity is 0
                  {
                    //if (in.position[i][1]<10e3) cout << "below S" <<"\n";
                    if (c == peridotite_idx)
                      reaction_rate_out->reaction_rates[i][c] =
                          old_porosity[i] * (old_peridotiteF[i] - old_peridotite[i]) / melting_time_scale;
                    else if (c == peridotiteF_idx)
                      reaction_rate_out->reaction_rates[i][c] =
                          old_porosity[i] * (old_peridotiteF[i] - old_peridotite[i]) / melting_time_scale; //
                    else if (c == porosity_idx)
                      reaction_rate_out->reaction_rates[i][c] =
                          -old_porosity[i] / melting_time_scale;
                    else
                      reaction_rate_out->reaction_rates[i][c] = 0.0;
                  }

                  else if (c_tot_old < c_f) // ... between solidus and liquidus
                  {
                    //if (in.position[i][1]<10e3) cout << "between SL" <<"\n";
                    const double Dc_s = c_s - old_peridotite[i];
                    const double Dc_f = c_f - old_peridotiteF[i];
                    // porosity update derived from lever rule
                    const double Dpor = ((1.0 - old_porosity[i]) * Dc_s + old_porosity[i] * Dc_f) /
                                        ((old_peridotite[i] + Dc_s / melting_time_scale * dtt) - (old_peridotiteF[i] + Dc_f / melting_time_scale * dtt));
                    if (c == peridotite_idx)
                      reaction_rate_out->reaction_rates[i][c] =
                          Dc_s / melting_time_scale;
                    else if (c == peridotiteF_idx)
                      reaction_rate_out->reaction_rates[i][c] =
                          Dc_f / melting_time_scale; // 
                    else if (c == porosity_idx)
                      reaction_rate_out->reaction_rates[i][c] = Dpor / melting_time_scale;
                    else
                      reaction_rate_out->reaction_rates[i][c] = 0.0;
                  }
                  else // if ( c_tot_old > c_f) // temperature (composition) above liquidus - equilibrium porosity is 1
                  {
                    //if (in.position[i][1]<10e3) cout << "above L" <<"\n";
                    if (c == peridotite_idx)
                      reaction_rate_out->reaction_rates[i][c] =
                          (1.0 - old_porosity[i]) * (old_peridotite[i] - old_peridotiteF[i]) / melting_time_scale;
                    else if (c == peridotiteF_idx)
                      reaction_rate_out->reaction_rates[i][c] =
                          (1.0 - old_porosity[i]) * (old_peridotite[i] - old_peridotiteF[i]) / melting_time_scale; // 
                    else if (c == porosity_idx)
                      reaction_rate_out->reaction_rates[i][c] =
                          (1.0 - old_porosity[i]) / melting_time_scale;
                    else
                      reaction_rate_out->reaction_rates[i][c] = 0.0;
                  }
                }
                out.reaction_terms[i][c] = 0.0;
              }
            }
          }
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
    MeltGlobalMod<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt global mod");
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
          prm.declare_entry("Reference melt viscosity", "10.",
                            Patterns::Double(0.),
                            "The value of the constant melt viscosity $\\eta_f$. Units: \\si{\\pascal\\second}.");
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
          prm.declare_entry("Composition density change", "0.0",
                            Patterns::Double(),
                            "Chenge of density due to compositional change of 1. "
                            "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
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
          prm.declare_entry("Exponential compositional strengthening factor", "0.0",
                            Patterns::Double(0.),
                            "Exponential dependency of solid viscosity on the (scaled) "
                            "solid composition (peridotite field). "
                            "Dimensionless factor. With a value of 0.0 (the default) the "
                            "viscosity does not depend on the composition.");
          prm.declare_entry("Maximum compositional viscosity change", "1e3",
                            Patterns::Double(0.),
                            "Maximum relative increase of solid viscosity due to composition variations. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MeltGlobalMod<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt global mod");
        {
          reference_rho_s = prm.get_double("Reference solid density");
          reference_rho_f = prm.get_double("Reference melt density");
          reference_T = prm.get_double("Reference temperature");
          eta_0 = prm.get_double("Reference shear viscosity");
          epsdot_0 = prm.get_double("Reference strainrate");
          eta_f = prm.get_double("Reference melt viscosity");
          reference_permeability = prm.get_double("Reference permeability");
          activation_energy = prm.get_double("Activation energy");
          stress_exponent = prm.get_double("Stress viscosity exponent");
          thermal_conductivity = prm.get_double("Thermal conductivity");
          reference_specific_heat = prm.get_double("Reference specific heat");
          thermal_expansivity = prm.get_double("Thermal expansion coefficient");
          alpha_phi = prm.get_double("Exponential melt weakening factor");
          composition_density_change = prm.get_double("Composition density change");
          surface_solidus = prm.get_double("Surface solidus");
          dT_solidus_liquidus = prm.get_double("Solidus liquidus difference");
          dC_solidus_liquidus = prm.get_double("Solidus liquidus composition difference");
          C_reference = prm.get_double("Reference solid composition");
          pressure_solidus_change = prm.get_double("Pressure solidus change");
          include_melting_and_freezing = prm.get_bool("Include melting and freezing");
          melting_time_scale = prm.get_double("Melting time scale for operator splitting");
          alpha_composition = prm.get_double("Exponential compositional strengthening factor");
          delta_eta_composition_max = prm.get_double("Maximum compositional viscosity change");

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
                        ExcMessage("Material model Melt global mod with melt transport only "
                                   "works if there is a compositional field called porosity."));
          }

          //if (this->include_melt_transport())
          {
            AssertThrow(this->get_parameters().use_operator_splitting,
                        ExcMessage("Material model Melt global mod with melt transport only "
                                   "works with operator splitting."));
          }
          AssertThrow(this->introspection().compositional_name_exists("porosity"),
                      ExcMessage("Material model Melt global mod only works if there is a "
                                 "compositional field called porosity."));
          AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                      ExcMessage("Material model Melt global mod only works if there is a "
                                 "compositional field called peridotite."));
          AssertThrow(this->introspection().compositional_name_exists("peridotiteF"),
                      ExcMessage("Material model Melt global mod only works if there is a "
                                 "compositional field called peridotiteF."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MeltGlobalMod<dim>::create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
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
    ASPECT_REGISTER_MATERIAL_MODEL(MeltGlobalMod,
                                   "melt global mod",
                                   "A material model for the modeling of melt transport. "
                                   "It includes a simplified melting parametrization with "
                                   "linear solidus and liquidus temperatures, "
                                   "and equilibrium melt fraction calculated by the lever rule.")
  }
}
