using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Xml.Linq;
using MedlineExtractor.Models;

namespace MedlineExtractor
{
    class Program
    {
        // Input and Output folder paths
        private static string _inputFolder = "";
        private static string _outputFolder = "";
        
        // Start and End of the file sequence
        private const int Start = 1;
        private const int End = 1219;
        
        // File name prefix
        private const string Prefix = "pubmed24n";
        
        // Headings to extract (in-order)
        // private static readonly string[] Headings = {"Id", "Title", "Abstract", "Country", "JournalName", "Year", "Mesh", "Authors"};
        private static readonly string[] Headings = {"PMID DOI", "Authors", "Journal.volume:page", "year", "title", "Keywords", "Abstract"};
        static void Main(string[] args)
        {
            if (args.Length > 0)
            {
                _inputFolder = args[0] ?? "";
                _outputFolder = args[1] ?? "";
            }
            
            for (var i = Start; i <= End; i++)
            {
                var fileNo = i.ToString("0000");
                var fileName = $"{Prefix}{fileNo}";
                var output = Process(fileName);
                
                Write(output, fileName);
                Console.WriteLine($"Completed {fileName}");
            }
        }

        private static IEnumerable<Output> Process(string fileName)
        {
            var fullPath = $"{_inputFolder}{fileName}.xml";
            var doc = XDocument.Load(fullPath);
            
            var articleElements = doc.Descendants(nameof(PubmedArticle));

            var outputs = new List<Output>();
            foreach (var record in articleElements)
            {
                var medlineCitation = record.Element(nameof(MedlineCitation));
                if (medlineCitation == null) continue;

                var article = medlineCitation.Element(nameof(Article));
                var authorList = medlineCitation
                    .Element(nameof(Article))
                    ?.Element(nameof(Article.AuthorList))?.Descendants(nameof(Author)).Select(author =>
                        new Author
                        {
                            ForeName = author.Element(nameof(Author.ForeName))?.Value,
                            LastName = author.Element(nameof(Author.LastName))?.Value
                        }).ToList();

                var keywordList = medlineCitation
                    ?.Element(nameof(MedlineCitation.KeywordList))?.Descendants(nameof(Keyword)).Select(keyword => keyword.Value).ToList();

                var meshHeadingList = medlineCitation
                    .Element(nameof(MedlineCitation.MeshHeadingList))
                    ?.Descendants(nameof(MeshHeading)).Select(mesh => new MeshHeading
                    {
                        DescriptorName = mesh.Element(nameof(MeshHeading.DescriptorName))?.Value
                    }).ToList();

                var output = new Output
                {
                    PubMedId = medlineCitation.Element(nameof(MedlineCitation.PMID))?.Value,
                    
                    Authors = string.Join(";",
                        (authorList ?? new List<Author>()).Select(x => x.ToString())),
                    
                    Keywords = string.Join(";", keywordList ?? new List<string>()),
                    
                    Mesh = "",

                    Year = article?.Element(nameof(Journal))
                        ?.Element(nameof(Journal.JournalIssue))
                        ?.Element(nameof(Journal.JournalIssue.PubDate))
                        ?.Element(nameof(Journal.JournalIssue.PubDate.Year))?.Value,
                    
                    Title = article?.Element(nameof(Article.ArticleTitle))?.Value,
                    
                    Abstract = article?.Element(nameof(Abstract))?.Element(nameof(Abstract.AbstractText))?.Value,
                    
                };
                outputs.Add(output);
            }

            return outputs;
        }

        private static void Write(IEnumerable<Output> output, string fileName)
        {
            using var outputFile = new StreamWriter(Path.Combine(_outputFolder, $"{fileName}.tsv"));
            outputFile.WriteLine(string.Join("\t", Headings));
            foreach (var line in output)
                outputFile.WriteLine(line);
        }
    }
}
