using System;
using System.Collections.Generic;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class ObsScalar
    {
        public int Id { get; set; }
        public int SimulationId { get; set; }
        public string ObsJson { get; set; }
        public string StatsObs { get; set; }
        public string NDat { get; set; }
        public string KDat { get; set; }
        public string SignDat { get; set; }
        public string DoccDat { get; set; }
        public DateTime? TsCreated { get; set; }

        public virtual Simulation Simulation { get; set; }
    }
}
