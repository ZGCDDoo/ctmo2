
[PSCredential]$credentials = Get-Credential -Message "Enter Legacy PostgreSQL database credentials for scaffolding"

[string] $userName = $credentials.UserName
[string] $passWord = ConvertFrom-SecureString -SecureString $credentials.Password -AsPlainText

[string] $provider = "Npgsql.EntityFrameworkCore.PostgreSQL"
[string] $connectionString = "Host=localhost;Database=ctmo2;Username=postgres;Password=$passWord"
Write-Host $connectionString


dotnet ef dbcontext scaffold $connectionString $provider              `
        --context SimulationLegacyContext --verbose                   `
        --startup-project "src/NCtmo2.API/NCtmo2.API.csproj"          `
        --project "src/NCtmo2.Core/NCtmo2.Core.csproj"                `
        --schema public --output-dir "src/NCtmo2.Core/Contexts/Legacy" 

