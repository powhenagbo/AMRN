# AMRN Deployment

## Frontend — Vercel
Use `amr-hgt-frontend` as the Vercel Root Directory.

- Framework: Vite
- Build command: `npm run build`
- Output directory: `dist`
- Environment variable: `VITE_API_BASE_URL=https://YOUR-BACKEND-URL`

The frontend falls back to `http://127.0.0.1:5050` during local development.

## Backend — Render
A starter `render.yaml` is included at the AMRN root. The Flask app now honors the platform `PORT` variable and `gunicorn` is included in `requirements.txt`.

External tools such as AMRFinderPlus/RGI may require additional system or Docker configuration on the backend host.
