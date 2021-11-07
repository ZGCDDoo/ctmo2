using System;
using System.Collections.Generic;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class SimulationsCompressed
    {
        public int Id { get; set; }
        public int SimulationId { get; set; }
        public byte[] Simulation { get; set; }

        public virtual Simulation SimulationNavigation { get; set; }
    }
}
