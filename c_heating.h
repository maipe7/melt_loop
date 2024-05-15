#ifndef _aspect_heating_model_composition_heating_h
#define _aspect_heating_model_composition_heating_h

#include <aspect/simulator_access.h>
#include <aspect/heating_model/interface.h>

namespace aspect
{
  namespace HeatingModel
  {

    template <int dim>
    class CompositionHeating : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /*
         * This heating model calculates heat production dependent on composition.
         * The composition is defined as
         * composition=(1-porosity)*peridotite+porosity*peridotiteF
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        double reference_heating;
        double reference_composition;
    };
  }
}

#endif