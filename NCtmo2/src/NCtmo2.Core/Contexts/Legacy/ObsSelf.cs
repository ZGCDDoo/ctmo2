using System;
using System.Collections.Generic;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class ObsSelf
    {
        public int Id { get; set; }
        public int SimulationId { get; set; }
        public string Conventions { get; set; }
        public byte[] SelfupMean { get; set; }
        public byte[] SelfupStd { get; set; }
        public byte[] SelfdownMean { get; set; }
        public byte[] SelfdownStd { get; set; }
        public DateTime? TsCreated { get; set; }

        public virtual Simulation Simulation { get; set; }
    }
}
