using System;
using System.Collections.Generic;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class ObsGf
    {
        public int Id { get; set; }
        public int SimulationId { get; set; }
        public string Conventions { get; set; }
        public byte[] GfupMean { get; set; }
        public byte[] GfupStd { get; set; }
        public byte[] GfdownMean { get; set; }
        public byte[] GfdownStd { get; set; }
        public DateTime? TsCreated { get; set; }

        public virtual Simulation Simulation { get; set; }
    }
}
