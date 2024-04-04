#include "c_heating.h"

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    CompositionHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {

        int ipor=this->introspection().compositional_index_for_name("porosity");
        int iper=this->introspection().compositional_index_for_name("peridotite");
        int iperF=this->introspection().compositional_index_for_name("peridotiteF");

        for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          double porosity = material_model_inputs.composition[q][ipor];
          double peridotite = material_model_inputs.composition[q][iper];
          double peridotiteF = material_model_inputs.composition[q][iperF];

          // The heating is H=H_{ref}*(1-(peridotite-porosity)) times density
          double composition = porosity*peridotiteF+(1.0-porosity)*peridotite;
          heating_model_outputs.heating_source_terms[q] = reference_heating
                                                       * composition/reference_composition ;
                                                       //* material_model_outputs.densities[q];
          heating_model_outputs.heating_source_terms[q] = std::max(heating_model_outputs.heating_source_terms[q],0.0);

          // the LHS heating term is zero
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }


    template <int dim>
    void
    CompositionHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Composition melt heating");
        {
          prm.declare_entry("Heat production for reference composition","0.0",
                            Patterns::Double(0.),
                            "Value of heat production for reference composition."
                            "Units: W/kg}.");
          prm.declare_entry("Reference composition","1.0",
                            Patterns::Double(0.),
                            "Reference composition."
                            "Units: --}.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    CompositionHeating<dim>::parse_parameters (ParameterHandler &prm)
    {

      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Composition melt heating");
        {
           reference_heating = prm.get_double ("Heat production for reference composition");
           reference_composition = prm.get_double ("Reference composition");
                     
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
      
      AssertThrow(reference_composition>0,
                          ExcMessage("Reference composition in Composition melt heating model must be positive!"));
      AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Porosity compositional field required in Composition melt heating model!"));
      AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                          ExcMessage("Peridotite compositional field required in Composition melt heating model!"));
      AssertThrow(this->introspection().compositional_name_exists("peridotiteF"),
                          ExcMessage("PeridotiteF compositional field required in Composition melt heating model!"));
      
    }


  }
}


namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(CompositionHeating,
                                  "composition melt heating",
                                  "Heating model where the heat production depends on"
                                  "composition. The composition is defined as:"
                                  "composition=(1-porosity)*peridotite+porosity*peridotiteF")
  }
}