
[PSCredential]$credentials = Get-Credential -Message "Enter Legacy PostgreSQL database credentials for scaffolding"

[string] $userName = $credentials.UserName
[string] $passWord = ConvertFrom-SecureString -SecureString $credentials.Password -AsPlainText

Write-Host $userName
Write-Host $passWord

