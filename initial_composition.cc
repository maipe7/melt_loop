
#include "initial_composition.h"

namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    double
    PPP<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const
    {
      AssertThrow(this->introspection().compositional_name_exists("porosity"),
                  ExcMessage("The initial composition plugin `PPP' did not find a "
                             "compositional field called `porosity' to initialize. Please add a "
                             "compositional field with this name."));
      AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                  ExcMessage("The initial composition plugin `PPP' did not find a "
                             "compositional field called `peridotite' to initialize. Please add a "
                             "compositional field with this name."));
      AssertThrow(this->introspection().compositional_name_exists("peridotiteF"),
                  ExcMessage("The initial composition plugin `PPP' did not find a "
                             "compositional field called `peridotiteF' to initialize. Please add a "
                             "compositional field with this name."));

      const unsigned int porosity_index = this->introspection().compositional_index_for_name("porosity");
      const unsigned int peridotite_index = this->introspection().compositional_index_for_name("peridotite");
      const unsigned int peridotiteF_index = this->introspection().compositional_index_for_name("peridotiteF");
      if (compositional_index == porosity_index)
      {
          // MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());

          // // Use the initial composition, except for the porosity, to prevent
          // // infinite recursion
          // for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
          //   if (i != porosity_index)
          //     in.composition[0][i] = this->get_initial_composition_manager().initial_composition(position,i);
          //   else
          //     in.composition[0][i] = 0.0;

          return 0.0; // zero porosity
      }
      else if (compositional_index ==  peridotite_index)
      {
        //double c0=0.7;
        double cvar0=0.1;
        double xsize=1e3;
        double ysize=10e3;
        double size0=200;
        int nanomaly=100;
        double x=position[0];
        double y=position[1];
        double c = c0;
        for (double i=0; i<nanomaly; i++)
        {
          srand(i);
          double xcenter=rand()%int(xsize);
          srand(i*2);
          double ycenter=rand()%int(ysize);

          srand(i*3);
          // allow both negative and positive anomalies
          double cvar=cvar0*0.1*(rand()%20-10);

          srand(i*4);
          double size=rand()%int(size0);

          c += cvar*std::exp(-(std::pow(x-xcenter,2)+std::pow(y-ycenter,2))/std::pow(size,2));
        }
        // composition diminishes to zero near bottom boundary
        if(y<1e3) c*= y/1e3;
        return c;
      }
      else if (compositional_index ==  peridotiteF_index)
      {
        // constant initial melt composition is valid for zero initial porosity
        return 2.0;
      }
      return 0.0;
    }

    template <int dim>
    void
    PPP<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("PPP");
        {
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
    PPP<dim>::parse_parameters (ParameterHandler &prm)
    {

      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("PPP");
        {
           c0 = prm.get_double ("Reference composition");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
      
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(PPP,
                                              "PPP",
                                              "TODO.")
  }
}
