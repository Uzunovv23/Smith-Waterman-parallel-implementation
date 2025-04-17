# Script for generating random sequences for testing
# Usage: .\generate_random_seq.ps1 -Length 10000 > ..\data\random_10k.fa

param (
    [Parameter(Mandatory=$true)][int]$Length,
    [string]$OutputFile = $null,
    [string]$Name = "Random sequence"
)

function Generate-RandomSequence {
    param (
        [int]$Length,
        [string]$Type = "DNA"
    )
    
    $bases = @('A', 'C', 'G', 'T')
    if ($Type -eq "Protein") {
        $bases = @('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
    }
    
    $random = New-Object System.Random
    $sequence = ""
    
    for ($i = 0; $i -lt $Length; $i++) {
        $base = $bases[$random.Next(0, $bases.Length)]
        $sequence += $base
        
        # New line every 80 characters for better readability
        if (($i+1) % 80 -eq 0) {
            $sequence += "`n"
        }
    }
    
    # Final new line if we didn't end on the 80th character
    if ($Length % 80 -ne 0 -and -not $sequence.EndsWith("`n")) {
        $sequence += "`n"
    }
    
    return $sequence
}

# Add a header in FASTA format
$output = ">$Name length $Length`n"

# Generate the sequence
$sequence = Generate-RandomSequence -Length $Length
$output += $sequence

# Output the result to a file or to the console
if ($OutputFile) {
    # Create the data directory if it doesn't exist
    $outputDir = Split-Path -Parent $OutputFile
    if ($outputDir -and -not (Test-Path $outputDir)) {
        New-Item -ItemType Directory -Path $outputDir -Force | Out-Null
    }
    
    $output | Out-File -FilePath $OutputFile -Encoding utf8
    Write-Host "The sequence has been generated and saved to: $OutputFile" -ForegroundColor Green
} else {
    # Output to console
    Write-Output $output
}

# Examples of usage:
# .\generate_random_seq.ps1 -Length 1000 > ..\data\random_1k.fa
# .\generate_random_seq.ps1 -Length 5000 -OutputFile ..\data\random_5k.fa -Name "Test DNA"