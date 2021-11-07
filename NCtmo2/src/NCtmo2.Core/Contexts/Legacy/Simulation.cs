using System;
using System.Collections.Generic;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class Simulation
    {
        public Simulation()
        {
            ConfigsBatches = new HashSet<ConfigsBatch>();
            ObsGfs = new HashSet<ObsGf>();
            ObsHybs = new HashSet<ObsHyb>();
            ObsScalars = new HashSet<ObsScalar>();
            ObsSelves = new HashSet<ObsSelf>();
            SimulationsCompresseds = new HashSet<SimulationsCompressed>();
        }

        public int Id { get; set; }
        public string ModelName { get; set; }
        public double U { get; set; }
        public double Mu { get; set; }
        public double Beta { get; set; }
        public double Uprime { get; set; }
        public double JH { get; set; }
        public int NOrb { get; set; }
        public double GPhonon { get; set; }
        public double W0Phonon { get; set; }
        public string CtmoVersion { get; set; }
        public string Params { get; set; }
        public DateTime? TsCreated { get; set; }

        public virtual ICollection<ConfigsBatch> ConfigsBatches { get; set; }
        public virtual ICollection<ObsGf> ObsGfs { get; set; }
        public virtual ICollection<ObsHyb> ObsHybs { get; set; }
        public virtual ICollection<ObsScalar> ObsScalars { get; set; }
        public virtual ICollection<ObsSelf> ObsSelves { get; set; }
        public virtual ICollection<SimulationsCompressed> SimulationsCompresseds { get; set; }
    }
}
