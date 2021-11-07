using System;
using System.Collections.Generic;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class TestObsHyb
    {
        public int Id { get; set; }
        public int SimulationId { get; set; }
        public string Conventions { get; set; }
        public byte[] HybupMean { get; set; }
        public byte[] HybupStd { get; set; }
        public byte[] HybdownMean { get; set; }
        public byte[] HybdownStd { get; set; }
        public DateTime? TsCreated { get; set; }

        public virtual TestSimulation Simulation { get; set; }
    }
}
