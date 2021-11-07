using System;
using Microsoft.EntityFrameworkCore;
using Microsoft.EntityFrameworkCore.Metadata;

#nullable disable

namespace NCtmo2.Core.Contexts.Legacy
{
    public partial class SimulationLegacyContext : DbContext
    {
        public SimulationLegacyContext()
        {
        }

        public SimulationLegacyContext(DbContextOptions<SimulationLegacyContext> options)
            : base(options)
        {
        }

        public virtual DbSet<ConfigsBatch> ConfigsBatches { get; set; }
        public virtual DbSet<ConfigsBatchesRaw> ConfigsBatchesRaws { get; set; }
        public virtual DbSet<ObsGf> ObsGfs { get; set; }
        public virtual DbSet<ObsHyb> ObsHybs { get; set; }
        public virtual DbSet<ObsScalar> ObsScalars { get; set; }
        public virtual DbSet<ObsSelf> ObsSelves { get; set; }
        public virtual DbSet<Simulation> Simulations { get; set; }
        public virtual DbSet<SimulationsCompressed> SimulationsCompresseds { get; set; }
        public virtual DbSet<TestConfigsBatch> TestConfigsBatches { get; set; }
        public virtual DbSet<TestObsGf> TestObsGfs { get; set; }
        public virtual DbSet<TestObsHyb> TestObsHybs { get; set; }
        public virtual DbSet<TestObsScalar> TestObsScalars { get; set; }
        public virtual DbSet<TestObsSelf> TestObsSelves { get; set; }
        public virtual DbSet<TestSimulation> TestSimulations { get; set; }
        public virtual DbSet<TestSimulationsCompressed> TestSimulationsCompresseds { get; set; }

        protected override void OnConfiguring(DbContextOptionsBuilder optionsBuilder)
        {
        }

        protected override void OnModelCreating(ModelBuilder modelBuilder)
        {
            modelBuilder.HasAnnotation("Relational:Collation", "English_United States.1252");

            modelBuilder.Entity<ConfigsBatch>(entity =>
            {
                entity.ToTable("configs__batches");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Batch)
                    .IsRequired()
                    .HasColumnName("batch");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.ConfigsBatches)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("configs__batches_simulation_id_fkey");
            });

            modelBuilder.Entity<ConfigsBatchesRaw>(entity =>
            {
                entity.ToTable("configs__batches__raw");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Batch)
                    .IsRequired()
                    .HasColumnName("batch");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");
            });

            modelBuilder.Entity<ObsGf>(entity =>
            {
                entity.ToTable("obs__gf");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Conventions).HasColumnName("conventions");

                entity.Property(e => e.GfdownMean).HasColumnName("gfdown_mean");

                entity.Property(e => e.GfdownStd).HasColumnName("gfdown_std");

                entity.Property(e => e.GfupMean)
                    .IsRequired()
                    .HasColumnName("gfup_mean");

                entity.Property(e => e.GfupStd)
                    .IsRequired()
                    .HasColumnName("gfup_std");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.ObsGfs)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("obs__gf_simulation_id_fkey");
            });

            modelBuilder.Entity<ObsHyb>(entity =>
            {
                entity.ToTable("obs__hyb");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Conventions).HasColumnName("conventions");

                entity.Property(e => e.HybdownMean).HasColumnName("hybdown_mean");

                entity.Property(e => e.HybdownStd).HasColumnName("hybdown_std");

                entity.Property(e => e.HybupMean)
                    .IsRequired()
                    .HasColumnName("hybup_mean");

                entity.Property(e => e.HybupStd)
                    .IsRequired()
                    .HasColumnName("hybup_std");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("now()");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.ObsHybs)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("obs__hyb_simulation_id_fkey");
            });

            modelBuilder.Entity<ObsScalar>(entity =>
            {
                entity.ToTable("obs__scalar");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.DoccDat)
                    .IsRequired()
                    .HasColumnName("docc_dat");

                entity.Property(e => e.KDat)
                    .IsRequired()
                    .HasColumnName("k_dat");

                entity.Property(e => e.NDat)
                    .IsRequired()
                    .HasColumnName("n_dat");

                entity.Property(e => e.ObsJson)
                    .IsRequired()
                    .HasColumnType("json")
                    .HasColumnName("obs_json");

                entity.Property(e => e.SignDat)
                    .IsRequired()
                    .HasColumnName("sign_dat");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.StatsObs)
                    .IsRequired()
                    .HasColumnType("json")
                    .HasColumnName("stats_obs");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.ObsScalars)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("obs__scalar_simulation_id_fkey");
            });

            modelBuilder.Entity<ObsSelf>(entity =>
            {
                entity.ToTable("obs__self");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Conventions).HasColumnName("conventions");

                entity.Property(e => e.SelfdownMean).HasColumnName("selfdown_mean");

                entity.Property(e => e.SelfdownStd).HasColumnName("selfdown_std");

                entity.Property(e => e.SelfupMean)
                    .IsRequired()
                    .HasColumnName("selfup_mean");

                entity.Property(e => e.SelfupStd)
                    .IsRequired()
                    .HasColumnName("selfup_std");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.ObsSelves)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("obs__self_simulation_id_fkey");
            });

            modelBuilder.Entity<Simulation>(entity =>
            {
                entity.ToTable("simulations");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Beta).HasColumnName("beta");

                entity.Property(e => e.CtmoVersion)
                    .IsRequired()
                    .HasMaxLength(128)
                    .HasColumnName("ctmo_version");

                entity.Property(e => e.GPhonon).HasColumnName("gPhonon");

                entity.Property(e => e.JH).HasColumnName("J_H");

                entity.Property(e => e.ModelName)
                    .HasMaxLength(128)
                    .HasColumnName("model_name");

                entity.Property(e => e.Mu).HasColumnName("mu");

                entity.Property(e => e.NOrb)
                    .HasColumnName("nOrb")
                    .HasDefaultValueSql("1");

                entity.Property(e => e.Params)
                    .IsRequired()
                    .HasColumnType("json")
                    .HasColumnName("params");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.Property(e => e.Uprime).HasColumnName("UPrime");

                entity.Property(e => e.W0Phonon).HasColumnName("w0Phonon");
            });

            modelBuilder.Entity<SimulationsCompressed>(entity =>
            {
                entity.ToTable("simulations__compressed");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Simulation)
                    .IsRequired()
                    .HasColumnName("simulation");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.HasOne(d => d.SimulationNavigation)
                    .WithMany(p => p.SimulationsCompresseds)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("simulations__compressed_simulation_id_fkey");
            });

            modelBuilder.Entity<TestConfigsBatch>(entity =>
            {
                entity.ToTable("test__configs__batches");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Batch)
                    .IsRequired()
                    .HasColumnName("batch");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.TestConfigsBatches)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("test__configs__batches_simulation_id_fkey");
            });

            modelBuilder.Entity<TestObsGf>(entity =>
            {
                entity.ToTable("test__obs__gf");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Conventions).HasColumnName("conventions");

                entity.Property(e => e.GfdownMean).HasColumnName("gfdown_mean");

                entity.Property(e => e.GfdownStd).HasColumnName("gfdown_std");

                entity.Property(e => e.GfupMean)
                    .IsRequired()
                    .HasColumnName("gfup_mean");

                entity.Property(e => e.GfupStd)
                    .IsRequired()
                    .HasColumnName("gfup_std");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.TestObsGfs)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("test__obs__gf_simulation_id_fkey");
            });

            modelBuilder.Entity<TestObsHyb>(entity =>
            {
                entity.ToTable("test__obs__hyb");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Conventions).HasColumnName("conventions");

                entity.Property(e => e.HybdownMean).HasColumnName("hybdown_mean");

                entity.Property(e => e.HybdownStd).HasColumnName("hybdown_std");

                entity.Property(e => e.HybupMean)
                    .IsRequired()
                    .HasColumnName("hybup_mean");

                entity.Property(e => e.HybupStd)
                    .IsRequired()
                    .HasColumnName("hybup_std");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("now()");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.TestObsHybs)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("test__obs__hyb_simulation_id_fkey");
            });

            modelBuilder.Entity<TestObsScalar>(entity =>
            {
                entity.ToTable("test__obs__scalar");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.DoccDat)
                    .IsRequired()
                    .HasColumnName("docc_dat");

                entity.Property(e => e.KDat)
                    .IsRequired()
                    .HasColumnName("k_dat");

                entity.Property(e => e.NDat)
                    .IsRequired()
                    .HasColumnName("n_dat");

                entity.Property(e => e.ObsJson)
                    .IsRequired()
                    .HasColumnType("json")
                    .HasColumnName("obs_json");

                entity.Property(e => e.SignDat)
                    .IsRequired()
                    .HasColumnName("sign_dat");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.StatsObs)
                    .IsRequired()
                    .HasColumnType("json")
                    .HasColumnName("stats_obs");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.TestObsScalars)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("test__obs__scalar_simulation_id_fkey");
            });

            modelBuilder.Entity<TestObsSelf>(entity =>
            {
                entity.ToTable("test__obs__self");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Conventions).HasColumnName("conventions");

                entity.Property(e => e.SelfdownMean).HasColumnName("selfdown_mean");

                entity.Property(e => e.SelfdownStd).HasColumnName("selfdown_std");

                entity.Property(e => e.SelfupMean)
                    .IsRequired()
                    .HasColumnName("selfup_mean");

                entity.Property(e => e.SelfupStd)
                    .IsRequired()
                    .HasColumnName("selfup_std");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.HasOne(d => d.Simulation)
                    .WithMany(p => p.TestObsSelves)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("test__obs__self_simulation_id_fkey");
            });

            modelBuilder.Entity<TestSimulation>(entity =>
            {
                entity.ToTable("test__simulations");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Beta).HasColumnName("beta");

                entity.Property(e => e.CtmoVersion)
                    .IsRequired()
                    .HasMaxLength(128)
                    .HasColumnName("ctmo_version");

                entity.Property(e => e.GPhonon).HasColumnName("gPhonon");

                entity.Property(e => e.JH).HasColumnName("J_H");

                entity.Property(e => e.ModelName)
                    .HasMaxLength(128)
                    .HasColumnName("model_name");

                entity.Property(e => e.Mu).HasColumnName("mu");

                entity.Property(e => e.NOrb)
                    .HasColumnName("nOrb")
                    .HasDefaultValueSql("1");

                entity.Property(e => e.Params)
                    .IsRequired()
                    .HasColumnType("json")
                    .HasColumnName("params");

                entity.Property(e => e.TsCreated)
                    .HasColumnType("timestamp with time zone")
                    .HasColumnName("ts_created")
                    .HasDefaultValueSql("CURRENT_TIMESTAMP");

                entity.Property(e => e.Uprime).HasColumnName("UPrime");

                entity.Property(e => e.W0Phonon).HasColumnName("w0Phonon");
            });

            modelBuilder.Entity<TestSimulationsCompressed>(entity =>
            {
                entity.ToTable("test__simulations__compressed");

                entity.Property(e => e.Id).HasColumnName("id");

                entity.Property(e => e.Simulation)
                    .IsRequired()
                    .HasColumnName("simulation");

                entity.Property(e => e.SimulationId).HasColumnName("simulation_id");

                entity.HasOne(d => d.SimulationNavigation)
                    .WithMany(p => p.TestSimulationsCompresseds)
                    .HasForeignKey(d => d.SimulationId)
                    .HasConstraintName("test__simulations__compressed_simulation_id_fkey");
            });

            OnModelCreatingPartial(modelBuilder);
        }

        partial void OnModelCreatingPartial(ModelBuilder modelBuilder);
    }
}
