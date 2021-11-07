using System;
using System.Collections.Generic;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class TestSimulation
    {
        public TestSimulation()
        {
            TestConfigsBatches = new HashSet<TestConfigsBatch>();
            TestObsGfs = new HashSet<TestObsGf>();
            TestObsHybs = new HashSet<TestObsHyb>();
            TestObsScalars = new HashSet<TestObsScalar>();
            TestObsSelves = new HashSet<TestObsSelf>();
            TestSimulationsCompresseds = new HashSet<TestSimulationsCompressed>();
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

        public virtual ICollection<TestConfigsBatch> TestConfigsBatches { get; set; }
        public virtual ICollection<TestObsGf> TestObsGfs { get; set; }
        public virtual ICollection<TestObsHyb> TestObsHybs { get; set; }
        public virtual ICollection<TestObsScalar> TestObsScalars { get; set; }
        public virtual ICollection<TestObsSelf> TestObsSelves { get; set; }
        public virtual ICollection<TestSimulationsCompressed> TestSimulationsCompresseds { get; set; }
    }
}
