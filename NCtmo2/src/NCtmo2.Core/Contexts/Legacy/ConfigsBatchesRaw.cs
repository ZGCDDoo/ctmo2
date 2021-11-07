using System;
using System.Collections.Generic;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class ConfigsBatchesRaw
    {
        public int Id { get; set; }
        public int SimulationId { get; set; }
        public byte[] Batch { get; set; }
        public DateTime? TsCreated { get; set; }
    }
}
