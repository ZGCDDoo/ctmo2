using System;
using System.Collections.Generic;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class ConfigsBatch
    {
        public int Id { get; set; }
        public int SimulationId { get; set; }
        public byte[] Batch { get; set; }
        public DateTime? TsCreated { get; set; }

        public virtual Simulation Simulation { get; set; }
    }
}
